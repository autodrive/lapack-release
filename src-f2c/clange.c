#line 1 "clange.f"
/* clange.f -- translated by f2c (version 20100827).
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

#line 1 "clange.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clange.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clange.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clange.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANGE( NORM, M, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, M, N */
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
/* > CLANGE  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANGE */
/* > \verbatim */
/* > */
/* >    CLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANGE as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0.  When M = 0, */
/* >          CLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0.  When N = 0, */
/* >          CLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The m by n matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
/* >          referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complexGEauxiliary */

/*  ===================================================================== */
doublereal clange_(char *norm, integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
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

#line 154 "clange.f"
    /* Parameter adjustments */
#line 154 "clange.f"
    a_dim1 = *lda;
#line 154 "clange.f"
    a_offset = 1 + a_dim1;
#line 154 "clange.f"
    a -= a_offset;
#line 154 "clange.f"
    --work;
#line 154 "clange.f"

#line 154 "clange.f"
    /* Function Body */
#line 154 "clange.f"
    if (min(*m,*n) == 0) {
#line 155 "clange.f"
	value = 0.;
#line 156 "clange.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 160 "clange.f"
	value = 0.;
#line 161 "clange.f"
	i__1 = *n;
#line 161 "clange.f"
	for (j = 1; j <= i__1; ++j) {
#line 162 "clange.f"
	    i__2 = *m;
#line 162 "clange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 163 "clange.f"
		temp = z_abs(&a[i__ + j * a_dim1]);
#line 164 "clange.f"
		if (value < temp || sisnan_(&temp)) {
#line 164 "clange.f"
		    value = temp;
#line 164 "clange.f"
		}
#line 165 "clange.f"
/* L10: */
#line 165 "clange.f"
	    }
#line 166 "clange.f"
/* L20: */
#line 166 "clange.f"
	}
#line 167 "clange.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 171 "clange.f"
	value = 0.;
#line 172 "clange.f"
	i__1 = *n;
#line 172 "clange.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "clange.f"
	    sum = 0.;
#line 174 "clange.f"
	    i__2 = *m;
#line 174 "clange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 175 "clange.f"
		sum += z_abs(&a[i__ + j * a_dim1]);
#line 176 "clange.f"
/* L30: */
#line 176 "clange.f"
	    }
#line 177 "clange.f"
	    if (value < sum || sisnan_(&sum)) {
#line 177 "clange.f"
		value = sum;
#line 177 "clange.f"
	    }
#line 178 "clange.f"
/* L40: */
#line 178 "clange.f"
	}
#line 179 "clange.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 183 "clange.f"
	i__1 = *m;
#line 183 "clange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 184 "clange.f"
	    work[i__] = 0.;
#line 185 "clange.f"
/* L50: */
#line 185 "clange.f"
	}
#line 186 "clange.f"
	i__1 = *n;
#line 186 "clange.f"
	for (j = 1; j <= i__1; ++j) {
#line 187 "clange.f"
	    i__2 = *m;
#line 187 "clange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 188 "clange.f"
		work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 189 "clange.f"
/* L60: */
#line 189 "clange.f"
	    }
#line 190 "clange.f"
/* L70: */
#line 190 "clange.f"
	}
#line 191 "clange.f"
	value = 0.;
#line 192 "clange.f"
	i__1 = *m;
#line 192 "clange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 193 "clange.f"
	    temp = work[i__];
#line 194 "clange.f"
	    if (value < temp || sisnan_(&temp)) {
#line 194 "clange.f"
		value = temp;
#line 194 "clange.f"
	    }
#line 195 "clange.f"
/* L80: */
#line 195 "clange.f"
	}
#line 196 "clange.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 200 "clange.f"
	scale = 0.;
#line 201 "clange.f"
	sum = 1.;
#line 202 "clange.f"
	i__1 = *n;
#line 202 "clange.f"
	for (j = 1; j <= i__1; ++j) {
#line 203 "clange.f"
	    classq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 204 "clange.f"
/* L90: */
#line 204 "clange.f"
	}
#line 205 "clange.f"
	value = scale * sqrt(sum);
#line 206 "clange.f"
    }

#line 208 "clange.f"
    ret_val = value;
#line 209 "clange.f"
    return ret_val;

/*     End of CLANGE */

} /* clange_ */

