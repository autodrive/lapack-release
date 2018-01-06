#line 1 "clanhs.f"
/* clanhs.f -- translated by f2c (version 20100827).
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

#line 1 "clanhs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of an upper Hessenberg matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLANHS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION CLANHS( NORM, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            LDA, N */
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
/* > CLANHS  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > Hessenberg matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANHS */
/* > \verbatim */
/* > */
/* >    CLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in CLANHS as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, CLANHS is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
doublereal clanhs_(char *norm, integer *n, doublecomplex *a, integer *lda, 
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
    extern /* Subroutine */ int classq_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    extern logical sisnan_(doublereal *);


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

#line 148 "clanhs.f"
    /* Parameter adjustments */
#line 148 "clanhs.f"
    a_dim1 = *lda;
#line 148 "clanhs.f"
    a_offset = 1 + a_dim1;
#line 148 "clanhs.f"
    a -= a_offset;
#line 148 "clanhs.f"
    --work;
#line 148 "clanhs.f"

#line 148 "clanhs.f"
    /* Function Body */
#line 148 "clanhs.f"
    if (*n == 0) {
#line 149 "clanhs.f"
	value = 0.;
#line 150 "clanhs.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 154 "clanhs.f"
	value = 0.;
#line 155 "clanhs.f"
	i__1 = *n;
#line 155 "clanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 156 "clanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 156 "clanhs.f"
	    i__2 = min(i__3,i__4);
#line 156 "clanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 157 "clanhs.f"
		sum = z_abs(&a[i__ + j * a_dim1]);
#line 158 "clanhs.f"
		if (value < sum || sisnan_(&sum)) {
#line 158 "clanhs.f"
		    value = sum;
#line 158 "clanhs.f"
		}
#line 159 "clanhs.f"
/* L10: */
#line 159 "clanhs.f"
	    }
#line 160 "clanhs.f"
/* L20: */
#line 160 "clanhs.f"
	}
#line 161 "clanhs.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 165 "clanhs.f"
	value = 0.;
#line 166 "clanhs.f"
	i__1 = *n;
#line 166 "clanhs.f"
	for (j = 1; j <= i__1; ++j) {
#line 167 "clanhs.f"
	    sum = 0.;
/* Computing MIN */
#line 168 "clanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 168 "clanhs.f"
	    i__2 = min(i__3,i__4);
#line 168 "clanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 169 "clanhs.f"
		sum += z_abs(&a[i__ + j * a_dim1]);
#line 170 "clanhs.f"
/* L30: */
#line 170 "clanhs.f"
	    }
#line 171 "clanhs.f"
	    if (value < sum || sisnan_(&sum)) {
#line 171 "clanhs.f"
		value = sum;
#line 171 "clanhs.f"
	    }
#line 172 "clanhs.f"
/* L40: */
#line 172 "clanhs.f"
	}
#line 173 "clanhs.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 177 "clanhs.f"
	i__1 = *n;
#line 177 "clanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "clanhs.f"
	    work[i__] = 0.;
#line 179 "clanhs.f"
/* L50: */
#line 179 "clanhs.f"
	}
#line 180 "clanhs.f"
	i__1 = *n;
#line 180 "clanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 181 "clanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 181 "clanhs.f"
	    i__2 = min(i__3,i__4);
#line 181 "clanhs.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 182 "clanhs.f"
		work[i__] += z_abs(&a[i__ + j * a_dim1]);
#line 183 "clanhs.f"
/* L60: */
#line 183 "clanhs.f"
	    }
#line 184 "clanhs.f"
/* L70: */
#line 184 "clanhs.f"
	}
#line 185 "clanhs.f"
	value = 0.;
#line 186 "clanhs.f"
	i__1 = *n;
#line 186 "clanhs.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 187 "clanhs.f"
	    sum = work[i__];
#line 188 "clanhs.f"
	    if (value < sum || sisnan_(&sum)) {
#line 188 "clanhs.f"
		value = sum;
#line 188 "clanhs.f"
	    }
#line 189 "clanhs.f"
/* L80: */
#line 189 "clanhs.f"
	}
#line 190 "clanhs.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 194 "clanhs.f"
	scale = 0.;
#line 195 "clanhs.f"
	sum = 1.;
#line 196 "clanhs.f"
	i__1 = *n;
#line 196 "clanhs.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 197 "clanhs.f"
	    i__3 = *n, i__4 = j + 1;
#line 197 "clanhs.f"
	    i__2 = min(i__3,i__4);
#line 197 "clanhs.f"
	    classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 198 "clanhs.f"
/* L90: */
#line 198 "clanhs.f"
	}
#line 199 "clanhs.f"
	value = scale * sqrt(sum);
#line 200 "clanhs.f"
    }

#line 202 "clanhs.f"
    ret_val = value;
#line 203 "clanhs.f"
    return ret_val;

/*     End of CLANHS */

} /* clanhs_ */

