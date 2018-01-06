#line 1 "slange.f"
/* slange.f -- translated by f2c (version 20100827).
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

#line 1 "slange.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slange.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slange.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slange.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANGE  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANGE */
/* > \verbatim */
/* > */
/* >    SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANGE as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0.  When M = 0, */
/* >          SLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0.  When N = 0, */
/* >          SLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \ingroup realGEauxiliary */

/*  ===================================================================== */
doublereal slange_(char *norm, integer *m, integer *n, doublereal *a, integer 
	*lda, doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
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

#line 152 "slange.f"
    /* Parameter adjustments */
#line 152 "slange.f"
    a_dim1 = *lda;
#line 152 "slange.f"
    a_offset = 1 + a_dim1;
#line 152 "slange.f"
    a -= a_offset;
#line 152 "slange.f"
    --work;
#line 152 "slange.f"

#line 152 "slange.f"
    /* Function Body */
#line 152 "slange.f"
    if (min(*m,*n) == 0) {
#line 153 "slange.f"
	value = 0.;
#line 154 "slange.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 158 "slange.f"
	value = 0.;
#line 159 "slange.f"
	i__1 = *n;
#line 159 "slange.f"
	for (j = 1; j <= i__1; ++j) {
#line 160 "slange.f"
	    i__2 = *m;
#line 160 "slange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 161 "slange.f"
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 162 "slange.f"
		if (value < temp || sisnan_(&temp)) {
#line 162 "slange.f"
		    value = temp;
#line 162 "slange.f"
		}
#line 163 "slange.f"
/* L10: */
#line 163 "slange.f"
	    }
#line 164 "slange.f"
/* L20: */
#line 164 "slange.f"
	}
#line 165 "slange.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 169 "slange.f"
	value = 0.;
#line 170 "slange.f"
	i__1 = *n;
#line 170 "slange.f"
	for (j = 1; j <= i__1; ++j) {
#line 171 "slange.f"
	    sum = 0.;
#line 172 "slange.f"
	    i__2 = *m;
#line 172 "slange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 173 "slange.f"
		sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 174 "slange.f"
/* L30: */
#line 174 "slange.f"
	    }
#line 175 "slange.f"
	    if (value < sum || sisnan_(&sum)) {
#line 175 "slange.f"
		value = sum;
#line 175 "slange.f"
	    }
#line 176 "slange.f"
/* L40: */
#line 176 "slange.f"
	}
#line 177 "slange.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 181 "slange.f"
	i__1 = *m;
#line 181 "slange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "slange.f"
	    work[i__] = 0.;
#line 183 "slange.f"
/* L50: */
#line 183 "slange.f"
	}
#line 184 "slange.f"
	i__1 = *n;
#line 184 "slange.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "slange.f"
	    i__2 = *m;
#line 185 "slange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 186 "slange.f"
		work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 187 "slange.f"
/* L60: */
#line 187 "slange.f"
	    }
#line 188 "slange.f"
/* L70: */
#line 188 "slange.f"
	}
#line 189 "slange.f"
	value = 0.;
#line 190 "slange.f"
	i__1 = *m;
#line 190 "slange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "slange.f"
	    temp = work[i__];
#line 192 "slange.f"
	    if (value < temp || sisnan_(&temp)) {
#line 192 "slange.f"
		value = temp;
#line 192 "slange.f"
	    }
#line 193 "slange.f"
/* L80: */
#line 193 "slange.f"
	}
#line 194 "slange.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 198 "slange.f"
	scale = 0.;
#line 199 "slange.f"
	sum = 1.;
#line 200 "slange.f"
	i__1 = *n;
#line 200 "slange.f"
	for (j = 1; j <= i__1; ++j) {
#line 201 "slange.f"
	    slassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 202 "slange.f"
/* L90: */
#line 202 "slange.f"
	}
#line 203 "slange.f"
	value = scale * sqrt(sum);
#line 204 "slange.f"
    }

#line 206 "slange.f"
    ret_val = value;
#line 207 "slange.f"
    return ret_val;

/*     End of SLANGE */

} /* slange_ */

