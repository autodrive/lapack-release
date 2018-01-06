#line 1 "dlange.f"
/* dlange.f -- translated by f2c (version 20100827).
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

#line 1 "dlange.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlange.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlange.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANGE  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANGE */
/* > \verbatim */
/* > */
/* >    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANGE as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0.  When M = 0, */
/* >          DLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0.  When N = 0, */
/* >          DLANGE is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
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

/* > \ingroup doubleGEauxiliary */

/*  ===================================================================== */
doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer 
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

#line 152 "dlange.f"
    /* Parameter adjustments */
#line 152 "dlange.f"
    a_dim1 = *lda;
#line 152 "dlange.f"
    a_offset = 1 + a_dim1;
#line 152 "dlange.f"
    a -= a_offset;
#line 152 "dlange.f"
    --work;
#line 152 "dlange.f"

#line 152 "dlange.f"
    /* Function Body */
#line 152 "dlange.f"
    if (min(*m,*n) == 0) {
#line 153 "dlange.f"
	value = 0.;
#line 154 "dlange.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 158 "dlange.f"
	value = 0.;
#line 159 "dlange.f"
	i__1 = *n;
#line 159 "dlange.f"
	for (j = 1; j <= i__1; ++j) {
#line 160 "dlange.f"
	    i__2 = *m;
#line 160 "dlange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 161 "dlange.f"
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 162 "dlange.f"
		if (value < temp || disnan_(&temp)) {
#line 162 "dlange.f"
		    value = temp;
#line 162 "dlange.f"
		}
#line 163 "dlange.f"
/* L10: */
#line 163 "dlange.f"
	    }
#line 164 "dlange.f"
/* L20: */
#line 164 "dlange.f"
	}
#line 165 "dlange.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 169 "dlange.f"
	value = 0.;
#line 170 "dlange.f"
	i__1 = *n;
#line 170 "dlange.f"
	for (j = 1; j <= i__1; ++j) {
#line 171 "dlange.f"
	    sum = 0.;
#line 172 "dlange.f"
	    i__2 = *m;
#line 172 "dlange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 173 "dlange.f"
		sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 174 "dlange.f"
/* L30: */
#line 174 "dlange.f"
	    }
#line 175 "dlange.f"
	    if (value < sum || disnan_(&sum)) {
#line 175 "dlange.f"
		value = sum;
#line 175 "dlange.f"
	    }
#line 176 "dlange.f"
/* L40: */
#line 176 "dlange.f"
	}
#line 177 "dlange.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 181 "dlange.f"
	i__1 = *m;
#line 181 "dlange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "dlange.f"
	    work[i__] = 0.;
#line 183 "dlange.f"
/* L50: */
#line 183 "dlange.f"
	}
#line 184 "dlange.f"
	i__1 = *n;
#line 184 "dlange.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "dlange.f"
	    i__2 = *m;
#line 185 "dlange.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 186 "dlange.f"
		work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 187 "dlange.f"
/* L60: */
#line 187 "dlange.f"
	    }
#line 188 "dlange.f"
/* L70: */
#line 188 "dlange.f"
	}
#line 189 "dlange.f"
	value = 0.;
#line 190 "dlange.f"
	i__1 = *m;
#line 190 "dlange.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "dlange.f"
	    temp = work[i__];
#line 192 "dlange.f"
	    if (value < temp || disnan_(&temp)) {
#line 192 "dlange.f"
		value = temp;
#line 192 "dlange.f"
	    }
#line 193 "dlange.f"
/* L80: */
#line 193 "dlange.f"
	}
#line 194 "dlange.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 198 "dlange.f"
	scale = 0.;
#line 199 "dlange.f"
	sum = 1.;
#line 200 "dlange.f"
	i__1 = *n;
#line 200 "dlange.f"
	for (j = 1; j <= i__1; ++j) {
#line 201 "dlange.f"
	    dlassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 202 "dlange.f"
/* L90: */
#line 202 "dlange.f"
	}
#line 203 "dlange.f"
	value = scale * sqrt(sum);
#line 204 "dlange.f"
    }

#line 206 "dlange.f"
    ret_val = value;
#line 207 "dlange.f"
    return ret_val;

/*     End of DLANGE */

} /* dlange_ */

