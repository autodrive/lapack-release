#line 1 "slangt.f"
/* slangt.f -- translated by f2c (version 20100827).
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

#line 1 "slangt.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slangt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slangt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slangt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANGT( NORM, N, DL, D, DU ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANGT  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANGT */
/* > \verbatim */
/* > */
/* >    SLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in SLANGT as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, SLANGT is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of A. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
doublereal slangt_(char *norm, integer *n, doublereal *dl, doublereal *d__, 
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

#line 144 "slangt.f"
    /* Parameter adjustments */
#line 144 "slangt.f"
    --du;
#line 144 "slangt.f"
    --d__;
#line 144 "slangt.f"
    --dl;
#line 144 "slangt.f"

#line 144 "slangt.f"
    /* Function Body */
#line 144 "slangt.f"
    if (*n <= 0) {
#line 145 "slangt.f"
	anorm = 0.;
#line 146 "slangt.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 150 "slangt.f"
	anorm = (d__1 = d__[*n], abs(d__1));
#line 151 "slangt.f"
	i__1 = *n - 1;
#line 151 "slangt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "slangt.f"
	    d__3 = (d__2 = dl[i__], abs(d__2));
#line 152 "slangt.f"
	    if (anorm < (d__1 = dl[i__], abs(d__1)) || sisnan_(&d__3)) {
#line 152 "slangt.f"
		anorm = (d__4 = dl[i__], abs(d__4));
#line 152 "slangt.f"
	    }
#line 154 "slangt.f"
	    d__3 = (d__2 = d__[i__], abs(d__2));
#line 154 "slangt.f"
	    if (anorm < (d__1 = d__[i__], abs(d__1)) || sisnan_(&d__3)) {
#line 154 "slangt.f"
		anorm = (d__4 = d__[i__], abs(d__4));
#line 154 "slangt.f"
	    }
#line 156 "slangt.f"
	    d__3 = (d__2 = du[i__], abs(d__2));
#line 156 "slangt.f"
	    if (anorm < (d__1 = du[i__], abs(d__1)) || sisnan_(&d__3)) {
#line 156 "slangt.f"
		anorm = (d__4 = du[i__], abs(d__4));
#line 156 "slangt.f"
	    }
#line 158 "slangt.f"
/* L10: */
#line 158 "slangt.f"
	}
#line 159 "slangt.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "slangt.f"
	if (*n == 1) {
#line 164 "slangt.f"
	    anorm = abs(d__[1]);
#line 165 "slangt.f"
	} else {
#line 166 "slangt.f"
	    anorm = abs(d__[1]) + abs(dl[1]);
#line 167 "slangt.f"
	    temp = (d__1 = d__[*n], abs(d__1)) + (d__2 = du[*n - 1], abs(d__2)
		    );
#line 168 "slangt.f"
	    if (anorm < temp || sisnan_(&temp)) {
#line 168 "slangt.f"
		anorm = temp;
#line 168 "slangt.f"
	    }
#line 169 "slangt.f"
	    i__1 = *n - 1;
#line 169 "slangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 170 "slangt.f"
		temp = (d__1 = d__[i__], abs(d__1)) + (d__2 = dl[i__], abs(
			d__2)) + (d__3 = du[i__ - 1], abs(d__3));
#line 171 "slangt.f"
		if (anorm < temp || sisnan_(&temp)) {
#line 171 "slangt.f"
		    anorm = temp;
#line 171 "slangt.f"
		}
#line 172 "slangt.f"
/* L20: */
#line 172 "slangt.f"
	    }
#line 173 "slangt.f"
	}
#line 174 "slangt.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 178 "slangt.f"
	if (*n == 1) {
#line 179 "slangt.f"
	    anorm = abs(d__[1]);
#line 180 "slangt.f"
	} else {
#line 181 "slangt.f"
	    anorm = abs(d__[1]) + abs(du[1]);
#line 182 "slangt.f"
	    temp = (d__1 = d__[*n], abs(d__1)) + (d__2 = dl[*n - 1], abs(d__2)
		    );
#line 183 "slangt.f"
	    if (anorm < temp || sisnan_(&temp)) {
#line 183 "slangt.f"
		anorm = temp;
#line 183 "slangt.f"
	    }
#line 184 "slangt.f"
	    i__1 = *n - 1;
#line 184 "slangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 185 "slangt.f"
		temp = (d__1 = d__[i__], abs(d__1)) + (d__2 = du[i__], abs(
			d__2)) + (d__3 = dl[i__ - 1], abs(d__3));
#line 186 "slangt.f"
		if (anorm < temp || sisnan_(&temp)) {
#line 186 "slangt.f"
		    anorm = temp;
#line 186 "slangt.f"
		}
#line 187 "slangt.f"
/* L30: */
#line 187 "slangt.f"
	    }
#line 188 "slangt.f"
	}
#line 189 "slangt.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 193 "slangt.f"
	scale = 0.;
#line 194 "slangt.f"
	sum = 1.;
#line 195 "slangt.f"
	slassq_(n, &d__[1], &c__1, &scale, &sum);
#line 196 "slangt.f"
	if (*n > 1) {
#line 197 "slangt.f"
	    i__1 = *n - 1;
#line 197 "slangt.f"
	    slassq_(&i__1, &dl[1], &c__1, &scale, &sum);
#line 198 "slangt.f"
	    i__1 = *n - 1;
#line 198 "slangt.f"
	    slassq_(&i__1, &du[1], &c__1, &scale, &sum);
#line 199 "slangt.f"
	}
#line 200 "slangt.f"
	anorm = scale * sqrt(sum);
#line 201 "slangt.f"
    }

#line 203 "slangt.f"
    ret_val = anorm;
#line 204 "slangt.f"
    return ret_val;

/*     End of SLANGT */

} /* slangt_ */

