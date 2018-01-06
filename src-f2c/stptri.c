#line 1 "stptri.f"
/* stptri.f -- translated by f2c (version 20100827).
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

#line 1 "stptri.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b STPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPTRI( UPLO, DIAG, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPTRI computes the inverse of a real upper or lower triangular */
/* > matrix A stored in packed format. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangular matrix A, stored */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n. */
/* >          See below for further details. */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same packed storage format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular */
/* >                matrix is singular and its inverse can not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A triangular matrix A can be transferred to packed storage using one */
/* >  of the following program segments: */
/* > */
/* >  UPLO = 'U':                      UPLO = 'L': */
/* > */
/* >        JC = 1                           JC = 1 */
/* >        DO 2 J = 1, N                    DO 2 J = 1, N */
/* >           DO 1 I = 1, J                    DO 1 I = J, N */
/* >              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J) */
/* >      1    CONTINUE                    1    CONTINUE */
/* >           JC = JC + J                      JC = JC + N - J + 1 */
/* >      2 CONTINUE                       2 CONTINUE */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int stptri_(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical upper;
    extern /* Subroutine */ int stpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer jclast;
    static logical nounit;


/*  -- LAPACK computational routine (version 3.7.0) -- */
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 155 "stptri.f"
    /* Parameter adjustments */
#line 155 "stptri.f"
    --ap;
#line 155 "stptri.f"

#line 155 "stptri.f"
    /* Function Body */
#line 155 "stptri.f"
    *info = 0;
#line 156 "stptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "stptri.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 158 "stptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 159 "stptri.f"
	*info = -1;
#line 160 "stptri.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 161 "stptri.f"
	*info = -2;
#line 162 "stptri.f"
    } else if (*n < 0) {
#line 163 "stptri.f"
	*info = -3;
#line 164 "stptri.f"
    }
#line 165 "stptri.f"
    if (*info != 0) {
#line 166 "stptri.f"
	i__1 = -(*info);
#line 166 "stptri.f"
	xerbla_("STPTRI", &i__1, (ftnlen)6);
#line 167 "stptri.f"
	return 0;
#line 168 "stptri.f"
    }

/*     Check for singularity if non-unit. */

#line 172 "stptri.f"
    if (nounit) {
#line 173 "stptri.f"
	if (upper) {
#line 174 "stptri.f"
	    jj = 0;
#line 175 "stptri.f"
	    i__1 = *n;
#line 175 "stptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 176 "stptri.f"
		jj += *info;
#line 177 "stptri.f"
		if (ap[jj] == 0.) {
#line 177 "stptri.f"
		    return 0;
#line 177 "stptri.f"
		}
#line 179 "stptri.f"
/* L10: */
#line 179 "stptri.f"
	    }
#line 180 "stptri.f"
	} else {
#line 181 "stptri.f"
	    jj = 1;
#line 182 "stptri.f"
	    i__1 = *n;
#line 182 "stptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 183 "stptri.f"
		if (ap[jj] == 0.) {
#line 183 "stptri.f"
		    return 0;
#line 183 "stptri.f"
		}
#line 185 "stptri.f"
		jj = jj + *n - *info + 1;
#line 186 "stptri.f"
/* L20: */
#line 186 "stptri.f"
	    }
#line 187 "stptri.f"
	}
#line 188 "stptri.f"
	*info = 0;
#line 189 "stptri.f"
    }

#line 191 "stptri.f"
    if (upper) {

/*        Compute inverse of upper triangular matrix. */

#line 195 "stptri.f"
	jc = 1;
#line 196 "stptri.f"
	i__1 = *n;
#line 196 "stptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 197 "stptri.f"
	    if (nounit) {
#line 198 "stptri.f"
		ap[jc + j - 1] = 1. / ap[jc + j - 1];
#line 199 "stptri.f"
		ajj = -ap[jc + j - 1];
#line 200 "stptri.f"
	    } else {
#line 201 "stptri.f"
		ajj = -1.;
#line 202 "stptri.f"
	    }

/*           Compute elements 1:j-1 of j-th column. */

#line 206 "stptri.f"
	    i__2 = j - 1;
#line 206 "stptri.f"
	    stpmv_("Upper", "No transpose", diag, &i__2, &ap[1], &ap[jc], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 208 "stptri.f"
	    i__2 = j - 1;
#line 208 "stptri.f"
	    sscal_(&i__2, &ajj, &ap[jc], &c__1);
#line 209 "stptri.f"
	    jc += j;
#line 210 "stptri.f"
/* L30: */
#line 210 "stptri.f"
	}

#line 212 "stptri.f"
    } else {

/*        Compute inverse of lower triangular matrix. */

#line 216 "stptri.f"
	jc = *n * (*n + 1) / 2;
#line 217 "stptri.f"
	for (j = *n; j >= 1; --j) {
#line 218 "stptri.f"
	    if (nounit) {
#line 219 "stptri.f"
		ap[jc] = 1. / ap[jc];
#line 220 "stptri.f"
		ajj = -ap[jc];
#line 221 "stptri.f"
	    } else {
#line 222 "stptri.f"
		ajj = -1.;
#line 223 "stptri.f"
	    }
#line 224 "stptri.f"
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

#line 228 "stptri.f"
		i__1 = *n - j;
#line 228 "stptri.f"
		stpmv_("Lower", "No transpose", diag, &i__1, &ap[jclast], &ap[
			jc + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 230 "stptri.f"
		i__1 = *n - j;
#line 230 "stptri.f"
		sscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
#line 231 "stptri.f"
	    }
#line 232 "stptri.f"
	    jclast = jc;
#line 233 "stptri.f"
	    jc = jc - *n + j - 2;
#line 234 "stptri.f"
/* L40: */
#line 234 "stptri.f"
	}
#line 235 "stptri.f"
    }

#line 237 "stptri.f"
    return 0;

/*     End of STPTRI */

} /* stptri_ */

