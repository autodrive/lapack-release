#line 1 "dlangt.f"
/* dlangt.f -- translated by f2c (version 20100827).
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

#line 1 "dlangt.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlangt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlangt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlangt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANGT( NORM, N, DL, D, DU ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANGT  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANGT */
/* > \verbatim */
/* > */
/* >    DLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANGT as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANGT is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of A. */
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
doublereal dlangt_(char *norm, integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum, temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
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

/*  ===================================================================== */

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

#line 144 "dlangt.f"
    /* Parameter adjustments */
#line 144 "dlangt.f"
    --du;
#line 144 "dlangt.f"
    --d__;
#line 144 "dlangt.f"
    --dl;
#line 144 "dlangt.f"

#line 144 "dlangt.f"
    /* Function Body */
#line 144 "dlangt.f"
    if (*n <= 0) {
#line 145 "dlangt.f"
	anorm = 0.;
#line 146 "dlangt.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 150 "dlangt.f"
	anorm = (d__1 = d__[*n], abs(d__1));
#line 151 "dlangt.f"
	i__1 = *n - 1;
#line 151 "dlangt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "dlangt.f"
	    d__3 = (d__2 = dl[i__], abs(d__2));
#line 152 "dlangt.f"
	    if (anorm < (d__1 = dl[i__], abs(d__1)) || disnan_(&d__3)) {
#line 152 "dlangt.f"
		anorm = (d__4 = dl[i__], abs(d__4));
#line 152 "dlangt.f"
	    }
#line 154 "dlangt.f"
	    d__3 = (d__2 = d__[i__], abs(d__2));
#line 154 "dlangt.f"
	    if (anorm < (d__1 = d__[i__], abs(d__1)) || disnan_(&d__3)) {
#line 154 "dlangt.f"
		anorm = (d__4 = d__[i__], abs(d__4));
#line 154 "dlangt.f"
	    }
#line 156 "dlangt.f"
	    d__3 = (d__2 = du[i__], abs(d__2));
#line 156 "dlangt.f"
	    if (anorm < (d__1 = du[i__], abs(d__1)) || disnan_(&d__3)) {
#line 156 "dlangt.f"
		anorm = (d__4 = du[i__], abs(d__4));
#line 156 "dlangt.f"
	    }
#line 158 "dlangt.f"
/* L10: */
#line 158 "dlangt.f"
	}
#line 159 "dlangt.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "dlangt.f"
	if (*n == 1) {
#line 164 "dlangt.f"
	    anorm = abs(d__[1]);
#line 165 "dlangt.f"
	} else {
#line 166 "dlangt.f"
	    anorm = abs(d__[1]) + abs(dl[1]);
#line 167 "dlangt.f"
	    temp = (d__1 = d__[*n], abs(d__1)) + (d__2 = du[*n - 1], abs(d__2)
		    );
#line 168 "dlangt.f"
	    if (anorm < temp || disnan_(&temp)) {
#line 168 "dlangt.f"
		anorm = temp;
#line 168 "dlangt.f"
	    }
#line 169 "dlangt.f"
	    i__1 = *n - 1;
#line 169 "dlangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 170 "dlangt.f"
		temp = (d__1 = d__[i__], abs(d__1)) + (d__2 = dl[i__], abs(
			d__2)) + (d__3 = du[i__ - 1], abs(d__3));
#line 171 "dlangt.f"
		if (anorm < temp || disnan_(&temp)) {
#line 171 "dlangt.f"
		    anorm = temp;
#line 171 "dlangt.f"
		}
#line 172 "dlangt.f"
/* L20: */
#line 172 "dlangt.f"
	    }
#line 173 "dlangt.f"
	}
#line 174 "dlangt.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 178 "dlangt.f"
	if (*n == 1) {
#line 179 "dlangt.f"
	    anorm = abs(d__[1]);
#line 180 "dlangt.f"
	} else {
#line 181 "dlangt.f"
	    anorm = abs(d__[1]) + abs(du[1]);
#line 182 "dlangt.f"
	    temp = (d__1 = d__[*n], abs(d__1)) + (d__2 = dl[*n - 1], abs(d__2)
		    );
#line 183 "dlangt.f"
	    if (anorm < temp || disnan_(&temp)) {
#line 183 "dlangt.f"
		anorm = temp;
#line 183 "dlangt.f"
	    }
#line 184 "dlangt.f"
	    i__1 = *n - 1;
#line 184 "dlangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 185 "dlangt.f"
		temp = (d__1 = d__[i__], abs(d__1)) + (d__2 = du[i__], abs(
			d__2)) + (d__3 = dl[i__ - 1], abs(d__3));
#line 186 "dlangt.f"
		if (anorm < temp || disnan_(&temp)) {
#line 186 "dlangt.f"
		    anorm = temp;
#line 186 "dlangt.f"
		}
#line 187 "dlangt.f"
/* L30: */
#line 187 "dlangt.f"
	    }
#line 188 "dlangt.f"
	}
#line 189 "dlangt.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 193 "dlangt.f"
	scale = 0.;
#line 194 "dlangt.f"
	sum = 1.;
#line 195 "dlangt.f"
	dlassq_(n, &d__[1], &c__1, &scale, &sum);
#line 196 "dlangt.f"
	if (*n > 1) {
#line 197 "dlangt.f"
	    i__1 = *n - 1;
#line 197 "dlangt.f"
	    dlassq_(&i__1, &dl[1], &c__1, &scale, &sum);
#line 198 "dlangt.f"
	    i__1 = *n - 1;
#line 198 "dlangt.f"
	    dlassq_(&i__1, &du[1], &c__1, &scale, &sum);
#line 199 "dlangt.f"
	}
#line 200 "dlangt.f"
	anorm = scale * sqrt(sum);
#line 201 "dlangt.f"
    }

#line 203 "dlangt.f"
    ret_val = anorm;
#line 204 "dlangt.f"
    return ret_val;

/*     End of DLANGT */

} /* dlangt_ */

