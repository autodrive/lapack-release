#line 1 "zlanht.f"
/* zlanht.f -- translated by f2c (version 20100827).
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

#line 1 "zlanht.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANHT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a complex Hermitian tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANHT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanht.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanht.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanht.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANHT( NORM, N, D, E ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ) */
/*       COMPLEX*16         E( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANHT  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex Hermitian tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANHT */
/* > \verbatim */
/* > */
/* >    ZLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANHT as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANHT is */
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
/* >          E is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) sub-diagonal or super-diagonal elements of A. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
doublereal zlanht_(char *norm, integer *n, doublereal *d__, doublecomplex *e, 
	ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *), zlassq_(integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 140 "zlanht.f"
    /* Parameter adjustments */
#line 140 "zlanht.f"
    --e;
#line 140 "zlanht.f"
    --d__;
#line 140 "zlanht.f"

#line 140 "zlanht.f"
    /* Function Body */
#line 140 "zlanht.f"
    if (*n <= 0) {
#line 141 "zlanht.f"
	anorm = 0.;
#line 142 "zlanht.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 146 "zlanht.f"
	anorm = (d__1 = d__[*n], abs(d__1));
#line 147 "zlanht.f"
	i__1 = *n - 1;
#line 147 "zlanht.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 148 "zlanht.f"
	    sum = (d__1 = d__[i__], abs(d__1));
#line 149 "zlanht.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 149 "zlanht.f"
		anorm = sum;
#line 149 "zlanht.f"
	    }
#line 150 "zlanht.f"
	    sum = z_abs(&e[i__]);
#line 151 "zlanht.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 151 "zlanht.f"
		anorm = sum;
#line 151 "zlanht.f"
	    }
#line 152 "zlanht.f"
/* L10: */
#line 152 "zlanht.f"
	}
#line 153 "zlanht.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find norm1(A). */

#line 158 "zlanht.f"
	if (*n == 1) {
#line 159 "zlanht.f"
	    anorm = abs(d__[1]);
#line 160 "zlanht.f"
	} else {
#line 161 "zlanht.f"
	    anorm = abs(d__[1]) + z_abs(&e[1]);
#line 162 "zlanht.f"
	    sum = z_abs(&e[*n - 1]) + (d__1 = d__[*n], abs(d__1));
#line 163 "zlanht.f"
	    if (anorm < sum || disnan_(&sum)) {
#line 163 "zlanht.f"
		anorm = sum;
#line 163 "zlanht.f"
	    }
#line 164 "zlanht.f"
	    i__1 = *n - 1;
#line 164 "zlanht.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 165 "zlanht.f"
		sum = (d__1 = d__[i__], abs(d__1)) + z_abs(&e[i__]) + z_abs(&
			e[i__ - 1]);
#line 166 "zlanht.f"
		if (anorm < sum || disnan_(&sum)) {
#line 166 "zlanht.f"
		    anorm = sum;
#line 166 "zlanht.f"
		}
#line 167 "zlanht.f"
/* L20: */
#line 167 "zlanht.f"
	    }
#line 168 "zlanht.f"
	}
#line 169 "zlanht.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 173 "zlanht.f"
	scale = 0.;
#line 174 "zlanht.f"
	sum = 1.;
#line 175 "zlanht.f"
	if (*n > 1) {
#line 176 "zlanht.f"
	    i__1 = *n - 1;
#line 176 "zlanht.f"
	    zlassq_(&i__1, &e[1], &c__1, &scale, &sum);
#line 177 "zlanht.f"
	    sum *= 2;
#line 178 "zlanht.f"
	}
#line 179 "zlanht.f"
	dlassq_(n, &d__[1], &c__1, &scale, &sum);
#line 180 "zlanht.f"
	anorm = scale * sqrt(sum);
#line 181 "zlanht.f"
    }

#line 183 "zlanht.f"
    ret_val = anorm;
#line 184 "zlanht.f"
    return ret_val;

/*     End of ZLANHT */

} /* zlanht_ */

