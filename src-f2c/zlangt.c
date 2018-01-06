#line 1 "zlangt.f"
/* zlangt.f -- translated by f2c (version 20100827).
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

#line 1 "zlangt.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute 
value of any element of a general tridiagonal matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlangt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlangt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlangt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION ZLANGT( NORM, N, DL, D, DU ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         D( * ), DL( * ), DU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANGT  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > complex tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return ZLANGT */
/* > \verbatim */
/* > */
/* >    ZLANGT = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* >          Specifies the value to be returned in ZLANGT as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, ZLANGT is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX*16 array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) super-diagonal elements of A. */
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
doublereal zlangt_(char *norm, integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sum, temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
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

#line 144 "zlangt.f"
    /* Parameter adjustments */
#line 144 "zlangt.f"
    --du;
#line 144 "zlangt.f"
    --d__;
#line 144 "zlangt.f"
    --dl;
#line 144 "zlangt.f"

#line 144 "zlangt.f"
    /* Function Body */
#line 144 "zlangt.f"
    if (*n <= 0) {
#line 145 "zlangt.f"
	anorm = 0.;
#line 146 "zlangt.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

#line 150 "zlangt.f"
	anorm = z_abs(&d__[*n]);
#line 151 "zlangt.f"
	i__1 = *n - 1;
#line 151 "zlangt.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "zlangt.f"
	    d__1 = z_abs(&dl[i__]);
#line 152 "zlangt.f"
	    if (anorm < z_abs(&dl[i__]) || disnan_(&d__1)) {
#line 152 "zlangt.f"
		anorm = z_abs(&dl[i__]);
#line 152 "zlangt.f"
	    }
#line 154 "zlangt.f"
	    d__1 = z_abs(&d__[i__]);
#line 154 "zlangt.f"
	    if (anorm < z_abs(&d__[i__]) || disnan_(&d__1)) {
#line 154 "zlangt.f"
		anorm = z_abs(&d__[i__]);
#line 154 "zlangt.f"
	    }
#line 156 "zlangt.f"
	    d__1 = z_abs(&du[i__]);
#line 156 "zlangt.f"
	    if (anorm < z_abs(&du[i__]) || disnan_(&d__1)) {
#line 156 "zlangt.f"
		anorm = z_abs(&du[i__]);
#line 156 "zlangt.f"
	    }
#line 158 "zlangt.f"
/* L10: */
#line 158 "zlangt.f"
	}
#line 159 "zlangt.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

#line 163 "zlangt.f"
	if (*n == 1) {
#line 164 "zlangt.f"
	    anorm = z_abs(&d__[1]);
#line 165 "zlangt.f"
	} else {
#line 166 "zlangt.f"
	    anorm = z_abs(&d__[1]) + z_abs(&dl[1]);
#line 167 "zlangt.f"
	    temp = z_abs(&d__[*n]) + z_abs(&du[*n - 1]);
#line 168 "zlangt.f"
	    if (anorm < temp || disnan_(&temp)) {
#line 168 "zlangt.f"
		anorm = temp;
#line 168 "zlangt.f"
	    }
#line 169 "zlangt.f"
	    i__1 = *n - 1;
#line 169 "zlangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 170 "zlangt.f"
		temp = z_abs(&d__[i__]) + z_abs(&dl[i__]) + z_abs(&du[i__ - 1]
			);
#line 171 "zlangt.f"
		if (anorm < temp || disnan_(&temp)) {
#line 171 "zlangt.f"
		    anorm = temp;
#line 171 "zlangt.f"
		}
#line 172 "zlangt.f"
/* L20: */
#line 172 "zlangt.f"
	    }
#line 173 "zlangt.f"
	}
#line 174 "zlangt.f"
    } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find normI(A). */

#line 178 "zlangt.f"
	if (*n == 1) {
#line 179 "zlangt.f"
	    anorm = z_abs(&d__[1]);
#line 180 "zlangt.f"
	} else {
#line 181 "zlangt.f"
	    anorm = z_abs(&d__[1]) + z_abs(&du[1]);
#line 182 "zlangt.f"
	    temp = z_abs(&d__[*n]) + z_abs(&dl[*n - 1]);
#line 183 "zlangt.f"
	    if (anorm < temp || disnan_(&temp)) {
#line 183 "zlangt.f"
		anorm = temp;
#line 183 "zlangt.f"
	    }
#line 184 "zlangt.f"
	    i__1 = *n - 1;
#line 184 "zlangt.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 185 "zlangt.f"
		temp = z_abs(&d__[i__]) + z_abs(&du[i__]) + z_abs(&dl[i__ - 1]
			);
#line 186 "zlangt.f"
		if (anorm < temp || disnan_(&temp)) {
#line 186 "zlangt.f"
		    anorm = temp;
#line 186 "zlangt.f"
		}
#line 187 "zlangt.f"
/* L30: */
#line 187 "zlangt.f"
	    }
#line 188 "zlangt.f"
	}
#line 189 "zlangt.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

#line 193 "zlangt.f"
	scale = 0.;
#line 194 "zlangt.f"
	sum = 1.;
#line 195 "zlangt.f"
	zlassq_(n, &d__[1], &c__1, &scale, &sum);
#line 196 "zlangt.f"
	if (*n > 1) {
#line 197 "zlangt.f"
	    i__1 = *n - 1;
#line 197 "zlangt.f"
	    zlassq_(&i__1, &dl[1], &c__1, &scale, &sum);
#line 198 "zlangt.f"
	    i__1 = *n - 1;
#line 198 "zlangt.f"
	    zlassq_(&i__1, &du[1], &c__1, &scale, &sum);
#line 199 "zlangt.f"
	}
#line 200 "zlangt.f"
	anorm = scale * sqrt(sum);
#line 201 "zlangt.f"
    }

#line 203 "zlangt.f"
    ret_val = anorm;
#line 204 "zlangt.f"
    return ret_val;

/*     End of ZLANGT */

} /* zlangt_ */

