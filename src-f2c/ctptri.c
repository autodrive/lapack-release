#line 1 "ctptri.f"
/* ctptri.f -- translated by f2c (version 20100827).
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

#line 1 "ctptri.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CTPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPTRI( UPLO, DIAG, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPTRI computes the inverse of a complex upper or lower triangular */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
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

/* > \ingroup complexOTHERcomputational */

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
/* Subroutine */ int ctptri_(char *uplo, char *diag, integer *n, 
	doublecomplex *ap, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, jc, jj;
    static doublecomplex ajj;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ctpmv_(char *, char *, char *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen, 
	    ftnlen);
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

#line 156 "ctptri.f"
    /* Parameter adjustments */
#line 156 "ctptri.f"
    --ap;
#line 156 "ctptri.f"

#line 156 "ctptri.f"
    /* Function Body */
#line 156 "ctptri.f"
    *info = 0;
#line 157 "ctptri.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 158 "ctptri.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 159 "ctptri.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 160 "ctptri.f"
	*info = -1;
#line 161 "ctptri.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 162 "ctptri.f"
	*info = -2;
#line 163 "ctptri.f"
    } else if (*n < 0) {
#line 164 "ctptri.f"
	*info = -3;
#line 165 "ctptri.f"
    }
#line 166 "ctptri.f"
    if (*info != 0) {
#line 167 "ctptri.f"
	i__1 = -(*info);
#line 167 "ctptri.f"
	xerbla_("CTPTRI", &i__1, (ftnlen)6);
#line 168 "ctptri.f"
	return 0;
#line 169 "ctptri.f"
    }

/*     Check for singularity if non-unit. */

#line 173 "ctptri.f"
    if (nounit) {
#line 174 "ctptri.f"
	if (upper) {
#line 175 "ctptri.f"
	    jj = 0;
#line 176 "ctptri.f"
	    i__1 = *n;
#line 176 "ctptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 177 "ctptri.f"
		jj += *info;
#line 178 "ctptri.f"
		i__2 = jj;
#line 178 "ctptri.f"
		if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
#line 178 "ctptri.f"
		    return 0;
#line 178 "ctptri.f"
		}
#line 180 "ctptri.f"
/* L10: */
#line 180 "ctptri.f"
	    }
#line 181 "ctptri.f"
	} else {
#line 182 "ctptri.f"
	    jj = 1;
#line 183 "ctptri.f"
	    i__1 = *n;
#line 183 "ctptri.f"
	    for (*info = 1; *info <= i__1; ++(*info)) {
#line 184 "ctptri.f"
		i__2 = jj;
#line 184 "ctptri.f"
		if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
#line 184 "ctptri.f"
		    return 0;
#line 184 "ctptri.f"
		}
#line 186 "ctptri.f"
		jj = jj + *n - *info + 1;
#line 187 "ctptri.f"
/* L20: */
#line 187 "ctptri.f"
	    }
#line 188 "ctptri.f"
	}
#line 189 "ctptri.f"
	*info = 0;
#line 190 "ctptri.f"
    }

#line 192 "ctptri.f"
    if (upper) {

/*        Compute inverse of upper triangular matrix. */

#line 196 "ctptri.f"
	jc = 1;
#line 197 "ctptri.f"
	i__1 = *n;
#line 197 "ctptri.f"
	for (j = 1; j <= i__1; ++j) {
#line 198 "ctptri.f"
	    if (nounit) {
#line 199 "ctptri.f"
		i__2 = jc + j - 1;
#line 199 "ctptri.f"
		z_div(&z__1, &c_b1, &ap[jc + j - 1]);
#line 199 "ctptri.f"
		ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
#line 200 "ctptri.f"
		i__2 = jc + j - 1;
#line 200 "ctptri.f"
		z__1.r = -ap[i__2].r, z__1.i = -ap[i__2].i;
#line 200 "ctptri.f"
		ajj.r = z__1.r, ajj.i = z__1.i;
#line 201 "ctptri.f"
	    } else {
#line 202 "ctptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 202 "ctptri.f"
		ajj.r = z__1.r, ajj.i = z__1.i;
#line 203 "ctptri.f"
	    }

/*           Compute elements 1:j-1 of j-th column. */

#line 207 "ctptri.f"
	    i__2 = j - 1;
#line 207 "ctptri.f"
	    ctpmv_("Upper", "No transpose", diag, &i__2, &ap[1], &ap[jc], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 209 "ctptri.f"
	    i__2 = j - 1;
#line 209 "ctptri.f"
	    cscal_(&i__2, &ajj, &ap[jc], &c__1);
#line 210 "ctptri.f"
	    jc += j;
#line 211 "ctptri.f"
/* L30: */
#line 211 "ctptri.f"
	}

#line 213 "ctptri.f"
    } else {

/*        Compute inverse of lower triangular matrix. */

#line 217 "ctptri.f"
	jc = *n * (*n + 1) / 2;
#line 218 "ctptri.f"
	for (j = *n; j >= 1; --j) {
#line 219 "ctptri.f"
	    if (nounit) {
#line 220 "ctptri.f"
		i__1 = jc;
#line 220 "ctptri.f"
		z_div(&z__1, &c_b1, &ap[jc]);
#line 220 "ctptri.f"
		ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
#line 221 "ctptri.f"
		i__1 = jc;
#line 221 "ctptri.f"
		z__1.r = -ap[i__1].r, z__1.i = -ap[i__1].i;
#line 221 "ctptri.f"
		ajj.r = z__1.r, ajj.i = z__1.i;
#line 222 "ctptri.f"
	    } else {
#line 223 "ctptri.f"
		z__1.r = -1., z__1.i = -0.;
#line 223 "ctptri.f"
		ajj.r = z__1.r, ajj.i = z__1.i;
#line 224 "ctptri.f"
	    }
#line 225 "ctptri.f"
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

#line 229 "ctptri.f"
		i__1 = *n - j;
#line 229 "ctptri.f"
		ctpmv_("Lower", "No transpose", diag, &i__1, &ap[jclast], &ap[
			jc + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
#line 231 "ctptri.f"
		i__1 = *n - j;
#line 231 "ctptri.f"
		cscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
#line 232 "ctptri.f"
	    }
#line 233 "ctptri.f"
	    jclast = jc;
#line 234 "ctptri.f"
	    jc = jc - *n + j - 2;
#line 235 "ctptri.f"
/* L40: */
#line 235 "ctptri.f"
	}
#line 236 "ctptri.f"
    }

#line 238 "ctptri.f"
    return 0;

/*     End of CTPTRI */

} /* ctptri_ */

