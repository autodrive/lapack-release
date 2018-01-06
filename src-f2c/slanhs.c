#line 1 "slanhs.f"
/* slanhs.f -- translated by f2c (version 20100827).
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

#line 1 "slanhs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of an upper Hessenberg matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANHS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slanhs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slanhs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slanhs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANHS( NORM, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANHS  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > Hessenberg matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANHS */
/* > \verbatim */
/* > */
/* >    SLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANHS as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
doublereal slanhs_(char *norm, integer *n, doublereal *a, integer *lda, 
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
    extern logical sisnan_(doublereal *);
    extern /* Subroutine */ int slassq_(integer *, doublereal *, integer *, 
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

#line 146 "slanhs.f"
    /* Parameter adjustments */
#line 146 "slanhs.f"
    a_dim1 = *lda;
#line 146 "slanhs.f"
    a_offset = 1 + a_dim1;
#line 146 "slanhs.f"
    a -= a_offset;
#line 146 "slanhs.f"
    --work;
#line 146 "slanhs.f"

#line 146 "slanhs.f"
    /* Function Body */
#line 146 "slanhs.f"
    if (*n == 0) {
#line 147 "slanhs.f"
	value = 0.;
#line 148 "slanhs.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 152 "slanhs.f"
	value = 0.;
#line 153 "slanhs.f"
	i__1 = *n;
#line 153 "slanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 154 "slanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 154 "slanhs.f"
	    i__2 = min(i__3,i__4);
#line 154 "slanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 155 "slanhs.f"
		sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 156 "slanhs.f"
		if (value < sum || sisnan_(&sum)) {
#line 156 "slanhs.f"
		    value = sum;
#line 156 "slanhs.f"
		}
#line 157 "slanhs.f"
/* L10: */
#line 157 "slanhs.f"
	    }
#line 158 "slanhs.f"
/* L20: */
#line 158 "slanhs.f"
	}
#line 159 "slanhs.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "slanhs.f"
	value = 0.;
#line 164 "slanhs.f"
	i__1 = *n;
#line 164 "slanhs.f"
	for (j = 1; j <= i__1; ++j) {
#line 165 "slanhs.f"
	    sum = 0.;
/* Computing MIN */
#line 166 "slanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 166 "slanhs.f"
	    i__2 = min(i__3,i__4);
#line 166 "slanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 167 "slanhs.f"
		sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 168 "slanhs.f"
/* L30: */
#line 168 "slanhs.f"
	    }
#line 169 "slanhs.f"
	    if (value < sum || sisnan_(&sum)) {
#line 169 "slanhs.f"
		value = sum;
#line 169 "slanhs.f"
	    }
#line 170 "slanhs.f"
/* L40: */
#line 170 "slanhs.f"
	}
#line 171 "slanhs.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 175 "slanhs.f"
	i__1 = *n;
#line 175 "slanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 176 "slanhs.f"
	    work[i__] = 0.;
#line 177 "slanhs.f"
/* L50: */
#line 177 "slanhs.f"
	}
#line 178 "slanhs.f"
	i__1 = *n;
#line 178 "slanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 179 "slanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 179 "slanhs.f"
	    i__2 = min(i__3,i__4);
#line 179 "slanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 180 "slanhs.f"
		work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 181 "slanhs.f"
/* L60: */
#line 181 "slanhs.f"
	    }
#line 182 "slanhs.f"
/* L70: */
#line 182 "slanhs.f"
	}
#line 183 "slanhs.f"
	value = 0.;
#line 184 "slanhs.f"
	i__1 = *n;
#line 184 "slanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "slanhs.f"
	    sum = work[i__];
#line 186 "slanhs.f"
	    if (value < sum || sisnan_(&sum)) {
#line 186 "slanhs.f"
		value = sum;
#line 186 "slanhs.f"
	    }
#line 187 "slanhs.f"
/* L80: */
#line 187 "slanhs.f"
	}
#line 188 "slanhs.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 192 "slanhs.f"
	scale = 0.;
#line 193 "slanhs.f"
	sum = 1.;
#line 194 "slanhs.f"
	i__1 = *n;
#line 194 "slanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 195 "slanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 195 "slanhs.f"
	    i__2 = min(i__3,i__4);
#line 195 "slanhs.f"
	    slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 196 "slanhs.f"
/* L90: */
#line 196 "slanhs.f"
	}
#line 197 "slanhs.f"
	value = scale * sqrt(sum);
#line 198 "slanhs.f"
    }

#line 200 "slanhs.f"
    ret_val = value;
#line 201 "slanhs.f"
    return ret_val;

/*     End of SLANHS */

} /* slanhs_ */

