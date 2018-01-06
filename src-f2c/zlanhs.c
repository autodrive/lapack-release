#line 1 "zlanhs.f"
/* zlanhs.f -- translated by f2c (version 20100827).
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

#line 1 "zlanhs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of an upper Hessenberg matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHS( NORM, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   WORK( * ) */
/*       COMPLEX*16         A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANHS  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > Hessenberg matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANHS */
/* > \verbatim */
/* > */
/* >    ZLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANHS as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANHS is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlanhs_(char *norm, integer *n, doublecomplex *a, integer *lda, 
	doublereal *work, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int zlassq_(integer *, doublecomplex *, integer *,
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 148 "zlanhs.f"
    /* Parameter adjustments */
#line 148 "zlanhs.f"
    a_dim1 = *lda;
#line 148 "zlanhs.f"
    a_offset = 1 + a_dim1;
#line 148 "zlanhs.f"
    a -= a_offset;
#line 148 "zlanhs.f"
    --work;
#line 148 "zlanhs.f"

#line 148 "zlanhs.f"
    /* Function Body */
#line 148 "zlanhs.f"
    if (*n == 0) {
#line 149 "zlanhs.f"
	value = 0.;
#line 150 "zlanhs.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 154 "zlanhs.f"
	value = 0.;
#line 155 "zlanhs.f"
	i__1 = *n;
#line 155 "zlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 156 "zlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 156 "zlanhs.f"
	    i__2 = min(i__3,i__4);
#line 156 "zlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 157 "zlanhs.f"
		sum = z_abs(&a[i__ + j * a_dim1]);
#line 158 "zlanhs.f"
		if (value < sum || disnan_(&sum)) {
#line 158 "zlanhs.f"
		    value = sum;
#line 158 "zlanhs.f"
		}
#line 159 "zlanhs.f"
/* L10: */
#line 159 "zlanhs.f"
	    }
#line 160 "zlanhs.f"
/* L20: */
#line 160 "zlanhs.f"
	}
#line 161 "zlanhs.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 165 "zlanhs.f"
	value = 0.;
#line 166 "zlanhs.f"
	i__1 = *n;
#line 166 "zlanhs.f"
	for (j = 1; j <= i__1; ++j) {
#line 167 "zlanhs.f"
	    sum = 0.;
/* Computing MIN */
#line 168 "zlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 168 "zlanhs.f"
	    i__2 = min(i__3,i__4);
#line 168 "zlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "zlanhs.f"
		sum += z_abs(&a[i__ + j * a_dim1]);
#line 170 "zlanhs.f"
/* L30: */
#line 170 "zlanhs.f"
	    }
#line 171 "zlanhs.f"
	    if (value < sum || disnan_(&sum)) {
#line 171 "zlanhs.f"
		value = sum;
#line 171 "zlanhs.f"
	    }
#line 172 "zlanhs.f"
/* L40: */
#line 172 "zlanhs.f"
	}
#line 173 "zlanhs.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 177 "zlanhs.f"
	i__1 = *n;
#line 177 "zlanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "zlanhs.f"
	    work[i__] = 0.;
#line 179 "zlanhs.f"
/* L50: */
#line 179 "zlanhs.f"
	}
#line 180 "zlanhs.f"
	i__1 = *n;
#line 180 "zlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 181 "zlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 181 "zlanhs.f"
	    i__2 = min(i__3,i__4);
#line 181 "zlanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 182 "zlanhs.f"
		work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 183 "zlanhs.f"
/* L60: */
#line 183 "zlanhs.f"
	    }
#line 184 "zlanhs.f"
/* L70: */
#line 184 "zlanhs.f"
	}
#line 185 "zlanhs.f"
	value = 0.;
#line 186 "zlanhs.f"
	i__1 = *n;
#line 186 "zlanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 187 "zlanhs.f"
	    sum = work[i__];
#line 188 "zlanhs.f"
	    if (value < sum || disnan_(&sum)) {
#line 188 "zlanhs.f"
		value = sum;
#line 188 "zlanhs.f"
	    }
#line 189 "zlanhs.f"
/* L80: */
#line 189 "zlanhs.f"
	}
#line 190 "zlanhs.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 194 "zlanhs.f"
	scale = 0.;
#line 195 "zlanhs.f"
	sum = 1.;
#line 196 "zlanhs.f"
	i__1 = *n;
#line 196 "zlanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 197 "zlanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 197 "zlanhs.f"
	    i__2 = min(i__3,i__4);
#line 197 "zlanhs.f"
	    zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 198 "zlanhs.f"
/* L90: */
#line 198 "zlanhs.f"
	}
#line 199 "zlanhs.f"
	value = scale * sqrt(sum);
#line 200 "zlanhs.f"
    }

#line 202 "zlanhs.f"
    ret_val = value;
#line 203 "zlanhs.f"
    return ret_val;

/*     End of ZLANHS */

} /* zlanhs_ */

