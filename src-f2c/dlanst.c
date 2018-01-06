#line 1 "dlanst.f"
/* dlanst.f -- translated by f2c (version 20100827).
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

#line 1 "dlanst.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a real symmetric tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANST  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANST */
/* > \verbatim */
/* > */
/* >    DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in DLANST as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANST is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) sub-diagonal or super-diagonal elements of A. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
doublereal dlanst_(char *norm, integer *n, doublereal *d__, doublereal *e, 
	ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum, scale;
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

#line 138 "dlanst.f"
    /* Parameter adjustments */
#line 138 "dlanst.f"
    --e;
#line 138 "dlanst.f"
    --d__;
#line 138 "dlanst.f"

#line 138 "dlanst.f"
    /* Function Body */
#line 138 "dlanst.f"
    if (*n <= 0) {
#line 139 "dlanst.f"
	anorm = 0.;
#line 140 "dlanst.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 144 "dlanst.f"
	anorm = (d__1 = d__[*n], abs(d__1));
#line 145 "dlanst.f"
	i__1 = *n - 1;
#line 145 "dlanst.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 146 "dlanst.f"
	    sum = (d__1 = d__[i__], abs(d__1));
#line 147 "dlanst.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 147 "dlanst.f"
		anorm = sum;
#line 147 "dlanst.f"
	    }
#line 148 "dlanst.f"
	    sum = (d__1 = e[i__], abs(d__1));
#line 149 "dlanst.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 149 "dlanst.f"
		anorm = sum;
#line 149 "dlanst.f"
	    }
#line 150 "dlanst.f"
/* L10: */
#line 150 "dlanst.f"
	}
#line 151 "dlanst.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find norm1(A). */

#line 156 "dlanst.f"
	if (*n == 1) {
#line 157 "dlanst.f"
	    anorm = abs(d__[1]);
#line 158 "dlanst.f"
	} else {
#line 159 "dlanst.f"
	    anorm = abs(d__[1]) + abs(e[1]);
#line 160 "dlanst.f"
	    sum = (d__1 = e[*n - 1], abs(d__1)) + (d__2 = d__[*n], abs(d__2));
#line 161 "dlanst.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 161 "dlanst.f"
		anorm = sum;
#line 161 "dlanst.f"
	    }
#line 162 "dlanst.f"
	    i__1 = *n - 1;
#line 162 "dlanst.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 163 "dlanst.f"
		sum = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[i__], abs(d__2)
			) + (d__3 = e[i__ - 1], abs(d__3));
#line 164 "dlanst.f"
		if (anorm < sum || disnan_(&sum)) {
#line 164 "dlanst.f"
		    anorm = sum;
#line 164 "dlanst.f"
		}
#line 165 "dlanst.f"
/* L20: */
#line 165 "dlanst.f"
	    }
#line 166 "dlanst.f"
	}
#line 167 "dlanst.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 171 "dlanst.f"
	scale = 0.;
#line 172 "dlanst.f"
	sum = 1.;
#line 173 "dlanst.f"
	if (*n > 1) {
#line 174 "dlanst.f"
	    i__1 = *n - 1;
#line 174 "dlanst.f"
	    dlassq_(&i__1, &e[1], &c__1, &scale, &sum);
#line 175 "dlanst.f"
	    sum *= 2;
#line 176 "dlanst.f"
	}
#line 177 "dlanst.f"
	dlassq_(n, &d__[1], &c__1, &scale, &sum);
#line 178 "dlanst.f"
	anorm = scale * sqrt(sum);
#line 179 "dlanst.f"
    }

#line 181 "dlanst.f"
    ret_val = anorm;
#line 182 "dlanst.f"
    return ret_val;

/*     End of DLANST */

} /* dlanst_ */

