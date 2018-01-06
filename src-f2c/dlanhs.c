#line 1 "dlanhs.f"
/* dlanhs.f -- translated by f2c (version 20100827).
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

#line 1 "dlanhs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of an upper Hessenberg matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANHS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanhs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanhs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanhs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
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
/* > DLANHS  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > Hessenberg matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANHS */
/* > \verbatim */
/* > */
/* >    DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANHS as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The n by n upper Hessenberg matrix A; the part of A below the */
/* >          first sub-diagonal is not referenced. */
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

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
doublereal dlanhs_(char *norm, integer *n, doublereal *a, integer *lda, 
	doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, scale;
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

#line 146 "dlanhs.f"
    /* Parameter adjustments */
#line 146 "dlanhs.f"
    a_dim1 = *lda;
#line 146 "dlanhs.f"
    a_offset = 1 + a_dim1;
#line 146 "dlanhs.f"
    a -= a_offset;
#line 146 "dlanhs.f"
    --work;
#line 146 "dlanhs.f"

#line 146 "dlanhs.f"
    /* Function Body */
#line 146 "dlanhs.f"
    if (*n == 0) {
#line 147 "dlanhs.f"
	value = 0.;
#line 148 "dlanhs.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 152 "dlanhs.f"
	value = 0.;
#line 153 "dlanhs.f"
	i__1 = *n;
#line 153 "dlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 154 "dlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 154 "dlanhs.f"
	    i__2 = min(i__3,i__4);
#line 154 "dlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 155 "dlanhs.f"
		sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 156 "dlanhs.f"
		if (value < sum || disnan_(&sum)) {
#line 156 "dlanhs.f"
		    value = sum;
#line 156 "dlanhs.f"
		}
#line 157 "dlanhs.f"
/* L10: */
#line 157 "dlanhs.f"
	    }
#line 158 "dlanhs.f"
/* L20: */
#line 158 "dlanhs.f"
	}
#line 159 "dlanhs.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "dlanhs.f"
	value = 0.;
#line 164 "dlanhs.f"
	i__1 = *n;
#line 164 "dlanhs.f"
	for (j = 1; j <= i__1; ++j) {
#line 165 "dlanhs.f"
	    sum = 0.;
/* Computing MIN */
#line 166 "dlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 166 "dlanhs.f"
	    i__2 = min(i__3,i__4);
#line 166 "dlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 167 "dlanhs.f"
		sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 168 "dlanhs.f"
/* L30: */
#line 168 "dlanhs.f"
	    }
#line 169 "dlanhs.f"
	    if (value < sum || disnan_(&sum)) {
#line 169 "dlanhs.f"
		value = sum;
#line 169 "dlanhs.f"
	    }
#line 170 "dlanhs.f"
/* L40: */
#line 170 "dlanhs.f"
	}
#line 171 "dlanhs.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 175 "dlanhs.f"
	i__1 = *n;
#line 175 "dlanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "dlanhs.f"
	    work[i__] = 0.;
#line 177 "dlanhs.f"
/* L50: */
#line 177 "dlanhs.f"
	}
#line 178 "dlanhs.f"
	i__1 = *n;
#line 178 "dlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 179 "dlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 179 "dlanhs.f"
	    i__2 = min(i__3,i__4);
#line 179 "dlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 180 "dlanhs.f"
		work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 181 "dlanhs.f"
/* L60: */
#line 181 "dlanhs.f"
	    }
#line 182 "dlanhs.f"
/* L70: */
#line 182 "dlanhs.f"
	}
#line 183 "dlanhs.f"
	value = 0.;
#line 184 "dlanhs.f"
	i__1 = *n;
#line 184 "dlanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "dlanhs.f"
	    sum = work[i__];
#line 186 "dlanhs.f"
	    if (value < sum || disnan_(&sum)) {
#line 186 "dlanhs.f"
		value = sum;
#line 186 "dlanhs.f"
	    }
#line 187 "dlanhs.f"
/* L80: */
#line 187 "dlanhs.f"
	}
#line 188 "dlanhs.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 192 "dlanhs.f"
	scale = 0.;
#line 193 "dlanhs.f"
	sum = 1.;
#line 194 "dlanhs.f"
	i__1 = *n;
#line 194 "dlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 195 "dlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 195 "dlanhs.f"
	    i__2 = min(i__3,i__4);
#line 195 "dlanhs.f"
	    dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 196 "dlanhs.f"
/* L90: */
#line 196 "dlanhs.f"
	}
#line 197 "dlanhs.f"
	value = scale * sqrt(sum);
#line 198 "dlanhs.f"
    }

#line 200 "dlanhs.f"
    ret_val = value;
#line 201 "dlanhs.f"
    return ret_val;

/*     End of DLANHS */

} /* dlanhs_ */

