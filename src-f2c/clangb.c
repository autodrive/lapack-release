#line 1 "clangb.f"
/* clangb.f -- translated by f2c (version 20100827).
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

#line 1 "clangb.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of general band matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clangb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clangb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clangb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANGB( NORM, N, KL, KU, AB, LDAB, */
/*                        WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            KL, KU, LDAB, N */
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
/* > CLANGB  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the element of  largest absolute value  of an */
/* > n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals. */
/* > \endverbatim */
/* > */
/* > \return CLANGB */
/* > \verbatim */
/* > */
/* >    CLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANGB as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANGB is */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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

/* > \ingroup complexGBauxiliary */

/*  ===================================================================== */
doublereal clangb_(char *norm, integer *n, integer *kl, integer *ku, 
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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 164 "clangb.f"
    /* Parameter adjustments */
#line 164 "clangb.f"
    ab_dim1 = *ldab;
#line 164 "clangb.f"
    ab_offset = 1 + ab_dim1;
#line 164 "clangb.f"
    ab -= ab_offset;
#line 164 "clangb.f"
    --work;
#line 164 "clangb.f"

#line 164 "clangb.f"
    /* Function Body */
#line 164 "clangb.f"
    if (*n == 0) {
#line 165 "clangb.f"
	value = 0.;
#line 166 "clangb.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 170 "clangb.f"
	value = 0.;
#line 171 "clangb.f"
	i__1 = *n;
#line 171 "clangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 172 "clangb.f"
	    i__2 = *ku + 2 - j;
/* Computing MIN */
#line 172 "clangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 172 "clangb.f"
	    i__3 = min(i__4,i__5);
#line 172 "clangb.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 173 "clangb.f"
		temp = z_abs(&ab[i__ + j * ab_dim1]);
#line 174 "clangb.f"
		if (value < temp || sisnan_(&temp)) {
#line 174 "clangb.f"
		    value = temp;
#line 174 "clangb.f"
		}
#line 175 "clangb.f"
/* L10: */
#line 175 "clangb.f"
	    }
#line 176 "clangb.f"
/* L20: */
#line 176 "clangb.f"
	}
#line 177 "clangb.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 181 "clangb.f"
	value = 0.;
#line 182 "clangb.f"
	i__1 = *n;
#line 182 "clangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 183 "clangb.f"
	    sum = 0.;
/* Computing MAX */
#line 184 "clangb.f"
	    i__3 = *ku + 2 - j;
/* Computing MIN */
#line 184 "clangb.f"
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 184 "clangb.f"
	    i__2 = min(i__4,i__5);
#line 184 "clangb.f"
	    for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
#line 185 "clangb.f"
		sum += z_abs(&ab[i__ + j * ab_dim1]);
#line 186 "clangb.f"
/* L30: */
#line 186 "clangb.f"
	    }
#line 187 "clangb.f"
	    if (value < sum || sisnan_(&sum)) {
#line 187 "clangb.f"
		value = sum;
#line 187 "clangb.f"
	    }
#line 188 "clangb.f"
/* L40: */
#line 188 "clangb.f"
	}
#line 189 "clangb.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 193 "clangb.f"
	i__1 = *n;
#line 193 "clangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "clangb.f"
	    work[i__] = 0.;
#line 195 "clangb.f"
/* L50: */
#line 195 "clangb.f"
	}
#line 196 "clangb.f"
	i__1 = *n;
#line 196 "clangb.f"
	for (j = 1; j <= i__1; ++j) {
#line 197 "clangb.f"
	    k = *ku + 1 - j;
/* Computing MAX */
#line 198 "clangb.f"
	    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
#line 198 "clangb.f"
	    i__5 = *n, i__6 = j + *kl;
#line 198 "clangb.f"
	    i__4 = min(i__5,i__6);
#line 198 "clangb.f"
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 199 "clangb.f"
		work[i__] += z_abs(&ab[k + i__ + j * ab_dim1]);
#line 200 "clangb.f"
/* L60: */
#line 200 "clangb.f"
	    }
#line 201 "clangb.f"
/* L70: */
#line 201 "clangb.f"
	}
#line 202 "clangb.f"
	value = 0.;
#line 203 "clangb.f"
	i__1 = *n;
#line 203 "clangb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 204 "clangb.f"
	    temp = work[i__];
#line 205 "clangb.f"
	    if (value < temp || sisnan_(&temp)) {
#line 205 "clangb.f"
		value = temp;
#line 205 "clangb.f"
	    }
#line 206 "clangb.f"
/* L80: */
#line 206 "clangb.f"
	}
#line 207 "clangb.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 211 "clangb.f"
	scale = 0.;
#line 212 "clangb.f"
	sum = 1.;
#line 213 "clangb.f"
	i__1 = *n;
#line 213 "clangb.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 214 "clangb.f"
	    i__4 = 1, i__2 = j - *ku;
#line 214 "clangb.f"
	    l = max(i__4,i__2);
#line 215 "clangb.f"
	    k = *ku + 1 - j + l;
/* Computing MIN */
#line 216 "clangb.f"
	    i__2 = *n, i__3 = j + *kl;
#line 216 "clangb.f"
	    i__4 = min(i__2,i__3) - l + 1;
#line 216 "clangb.f"
	    classq_(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
#line 217 "clangb.f"
/* L90: */
#line 217 "clangb.f"
	}
#line 218 "clangb.f"
	value = scale * sqrt(sum);
#line 219 "clangb.f"
    }

#line 221 "clangb.f"
    ret_val = value;
#line 222 "clangb.f"
    return ret_val;

/*     End of CLANGB */

} /* clangb_ */

