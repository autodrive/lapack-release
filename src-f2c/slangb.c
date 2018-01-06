#line 1 "slangb.f"
/* slangb.f -- translated by f2c (version 20100827).
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

#line 1 "slangb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of general band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slangb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slangb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slangb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANGB( NORM, N, KL, KU, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            KL, KU, LDAB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANGB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return SLANGB */
/* > \verbatim */
/* > */
/* >    SLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANGB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANGB is */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I'; otherwise, WORK is not */
/* >          referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGBauxiliary */

/*  ===================================================================== */
doublereal slangb_(char *norm, integer *n, integer *kl, integer *ku, 
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
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

#line 163 "slangb.f"
    /* Parameter adjustments */
#line 163 "slangb.f"
    ab_dim1 = *ldab;
#line 163 "slangb.f"
    ab_offset = 1 + ab_dim1;
#line 163 "slangb.f"
    ab -= ab_offset;
#line 163 "slangb.f"
    --work;
#line 163 "slangb.f"

#line 163 "slangb.f"
    /* Function Body */
#line 163 "slangb.f"
    if (*n == 0) {
#line 164 "slangb.f"
	value = 0.;
#line 165 "slangb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 169 "slangb.f"
	value = 0.;
#line 170 "slangb.f"
	i__1 = *n;
#line 170 "slangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 171 "slangb.f"
	    i__2 = *ku + 2 - j;
/* Computing MIN */
#line 171 "slangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 171 "slangb.f"
	    i__3 = min(i__4,i__5);
#line 171 "slangb.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 172 "slangb.f"
		temp = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 173 "slangb.f"
		if (value < temp || sisnan_(&temp)) {
#line 173 "slangb.f"
		    value = temp;
#line 173 "slangb.f"
		}
#line 174 "slangb.f"
/* L10: */
#line 174 "slangb.f"
	    }
#line 175 "slangb.f"
/* L20: */
#line 175 "slangb.f"
	}
#line 176 "slangb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 180 "slangb.f"
	value = 0.;
#line 181 "slangb.f"
	i__1 = *n;
#line 181 "slangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 182 "slangb.f"
	    sum = 0.;
/* Computing MAX */
#line 183 "slangb.f"
	    i__3 = *ku + 2 - j;
/* Computing MIN */
#line 183 "slangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 183 "slangb.f"
	    i__2 = min(i__4,i__5);
#line 183 "slangb.f"
	    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 184 "slangb.f"
		sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
#line 185 "slangb.f"
/* L30: */
#line 185 "slangb.f"
	    }
#line 186 "slangb.f"
	    if (value < sum || sisnan_(&sum)) {
#line 186 "slangb.f"
		value = sum;
#line 186 "slangb.f"
	    }
#line 187 "slangb.f"
/* L40: */
#line 187 "slangb.f"
	}
#line 188 "slangb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 192 "slangb.f"
	i__1 = *n;
#line 192 "slangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 193 "slangb.f"
	    work[i__] = 0.;
#line 194 "slangb.f"
/* L50: */
#line 194 "slangb.f"
	}
#line 195 "slangb.f"
	i__1 = *n;
#line 195 "slangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "slangb.f"
	    k = *ku + 1 - j;
/* Computing MAX */
#line 197 "slangb.f"
	    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 197 "slangb.f"
	    i__5 = *n, i__6 = j + *kl;
#line 197 "slangb.f"
	    i__4 = min(i__5,i__6);
#line 197 "slangb.f"
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 198 "slangb.f"
		work[i__] += (d__1 = ab[k + i__ + j * ab_dim1], abs(d__1));
#line 199 "slangb.f"
/* L60: */
#line 199 "slangb.f"
	    }
#line 200 "slangb.f"
/* L70: */
#line 200 "slangb.f"
	}
#line 201 "slangb.f"
	value = 0.;
#line 202 "slangb.f"
	i__1 = *n;
#line 202 "slangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "slangb.f"
	    temp = work[i__];
#line 204 "slangb.f"
	    if (value < temp || sisnan_(&temp)) {
#line 204 "slangb.f"
		value = temp;
#line 204 "slangb.f"
	    }
#line 205 "slangb.f"
/* L80: */
#line 205 "slangb.f"
	}
#line 206 "slangb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 210 "slangb.f"
	scale = 0.;
#line 211 "slangb.f"
	sum = 1.;
#line 212 "slangb.f"
	i__1 = *n;
#line 212 "slangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 213 "slangb.f"
	    i__4 = 1, i__2 = j - *ku;
#line 213 "slangb.f"
	    l = max(i__4,i__2);
#line 214 "slangb.f"
	    k = *ku + 1 - j + l;
/* Computing MIN */
#line 215 "slangb.f"
	    i__2 = *n, i__3 = j + *kl;
#line 215 "slangb.f"
	    i__4 = min(i__2,i__3) - l + 1;
#line 215 "slangb.f"
	    slassq_(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
#line 216 "slangb.f"
/* L90: */
#line 216 "slangb.f"
	}
#line 217 "slangb.f"
	value = scale * sqrt(sum);
#line 218 "slangb.f"
    }

#line 220 "slangb.f"
    ret_val = value;
#line 221 "slangb.f"
    return ret_val;

/*     End of SLANGB */

} /* slangb_ */

