#line 1 "dtptri.f"
/* dtptri.f -- translated by f2c (version 20100827).
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

#line 1 "dtptri.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DTPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPTRI computes the inverse of a real upper or lower triangular */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

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
/* Subroutine */ int dtptri_(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, jc, jj;
    static doublereal ajj;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer jclast;
    static logical nounit;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 155 "dtptri.f"
    /* Parameter adjustments */
#line 155 "dtptri.f"
    --ap;
#line 155 "dtptri.f"

#line 155 "dtptri.f"
    /* Function Body */
#line 155 "dtptri.f"
    *info = 0;
#line 156 "dtptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "dtptri.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 158 "dtptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 159 "dtptri.f"
	*info = -1;
#line 160 "dtptri.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 161 "dtptri.f"
	*info = -2;
#line 162 "dtptri.f"
    } else if (*n < 0) {
#line 163 "dtptri.f"
	*info = -3;
#line 164 "dtptri.f"
    }
#line 165 "dtptri.f"
    if (*info != 0) {
#line 166 "dtptri.f"
	i__1 = -(*info);
#line 166 "dtptri.f"
	xerbla_("DTPTRI", &i__1, (ftnlen)6);
#line 167 "dtptri.f"
	return 0;
#line 168 "dtptri.f"
    }

/*     Check for singularity if non-unit. */

#line 172 "dtptri.f"
    if (nounit) {
#line 173 "dtptri.f"
	if (upper) {
#line 174 "dtptri.f"
	    jj = 0;
#line 175 "dtptri.f"
	    i__1 = *n;
#line 175 "dtptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 176 "dtptri.f"
		jj += *info;
#line 177 "dtptri.f"
		if (ap[jj] == 0.) {
#line 177 "dtptri.f"
		    return 0;
#line 177 "dtptri.f"
		}
#line 179 "dtptri.f"
/* L10: */
#line 179 "dtptri.f"
	    }
#line 180 "dtptri.f"
	} else {
#line 181 "dtptri.f"
	    jj = 1;
#line 182 "dtptri.f"
	    i__1 = *n;
#line 182 "dtptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 183 "dtptri.f"
		if (ap[jj] == 0.) {
#line 183 "dtptri.f"
		    return 0;
#line 183 "dtptri.f"
		}
#line 185 "dtptri.f"
		jj = jj + *n - *info + 1;
#line 186 "dtptri.f"
/* L20: */
#line 186 "dtptri.f"
	    }
#line 187 "dtptri.f"
	}
#line 188 "dtptri.f"
	*info = 0;
#line 189 "dtptri.f"
    }

#line 191 "dtptri.f"
    if (upper) {

/*        Compute inverse of upper triangular matrix. */

#line 195 "dtptri.f"
	jc = 1;
#line 196 "dtptri.f"
	i__1 = *n;
#line 196 "dtptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 197 "dtptri.f"
	    if (nounit) {
#line 198 "dtptri.f"
		ap[jc + j - 1] = 1. / ap[jc + j - 1];
#line 199 "dtptri.f"
		ajj = -ap[jc + j - 1];
#line 200 "dtptri.f"
	    } else {
#line 201 "dtptri.f"
		ajj = -1.;
#line 202 "dtptri.f"
	    }

/*           Compute elements 1:j-1 of j-th column. */

#line 206 "dtptri.f"
	    i__2 = j - 1;
#line 206 "dtptri.f"
	    dtpmv_("Upper", "No transpose", diag, &i__2, &ap[1], &ap[jc], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 208 "dtptri.f"
	    i__2 = j - 1;
#line 208 "dtptri.f"
	    dscal_(&i__2, &ajj, &ap[jc], &c__1);
#line 209 "dtptri.f"
	    jc += j;
#line 210 "dtptri.f"
/* L30: */
#line 210 "dtptri.f"
	}

#line 212 "dtptri.f"
    } else {

/*        Compute inverse of lower triangular matrix. */

#line 216 "dtptri.f"
	jc = *n * (*n + 1) / 2;
#line 217 "dtptri.f"
	for (j = *n; j >= 1; --j) {
#line 218 "dtptri.f"
	    if (nounit) {
#line 219 "dtptri.f"
		ap[jc] = 1. / ap[jc];
#line 220 "dtptri.f"
		ajj = -ap[jc];
#line 221 "dtptri.f"
	    } else {
#line 222 "dtptri.f"
		ajj = -1.;
#line 223 "dtptri.f"
	    }
#line 224 "dtptri.f"
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

#line 228 "dtptri.f"
		i__1 = *n - j;
#line 228 "dtptri.f"
		dtpmv_("Lower", "No transpose", diag, &i__1, &ap[jclast], &ap[
			jc + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 230 "dtptri.f"
		i__1 = *n - j;
#line 230 "dtptri.f"
		dscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
#line 231 "dtptri.f"
	    }
#line 232 "dtptri.f"
	    jclast = jc;
#line 233 "dtptri.f"
	    jc = jc - *n + j - 2;
#line 234 "dtptri.f"
/* L40: */
#line 234 "dtptri.f"
	}
#line 235 "dtptri.f"
    }

#line 237 "dtptri.f"
    return 0;

/*     End of DTPTRI */

} /* dtptri_ */

