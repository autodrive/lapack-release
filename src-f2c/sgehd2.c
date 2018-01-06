#line 1 "sgehd2.f"
/* sgehd2.f -- translated by f2c (version 20100827).
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

#line 1 "sgehd2.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm. 
*/

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEHD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgehd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgehd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgehd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEHD2 reduces a real general matrix A to upper Hessenberg form H by */
/* > an orthogonal similarity transformation:  Q**T * A * Q = H . */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* > */
/* >          It is assumed that A is already upper triangular in rows */
/* >          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/* >          set by a previous call to SGEBAL; otherwise they should be */
/* >          set to 1 and N respectively. See Further Details. */
/* >          1 <= ILO <= IHI <= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the n by n general matrix to be reduced. */
/* >          On exit, the upper triangle and the first subdiagonal of A */
/* >          are overwritten with the upper Hessenberg matrix H, and the */
/* >          elements below the first subdiagonal, with the array TAU, */
/* >          represent the orthogonal matrix Q as a product of elementary */
/* >          reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix Q is represented as a product of (ihi-ilo) elementary */
/* >  reflectors */
/* > */
/* >     Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**T */
/* > */
/* >  where tau is a real scalar, and v is a real vector with */
/* >  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on */
/* >  exit in A(i+2:ihi,i), and tau in TAU(i). */
/* > */
/* >  The contents of A are illustrated by the following example, with */
/* >  n = 7, ilo = 2 and ihi = 6: */
/* > */
/* >  on entry,                        on exit, */
/* > */
/* >  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a ) */
/* >  (     a   a   a   a   a   a )    (      a   h   h   h   h   a ) */
/* >  (     a   a   a   a   a   a )    (      h   h   h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h ) */
/* >  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h ) */
/* >  (                         a )    (                          a ) */
/* > */
/* >  where a denotes an element of the original matrix A, h denotes a */
/* >  modified element of the upper Hessenberg matrix H, and vi denotes an */
/* >  element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sgehd2_(integer *n, integer *ilo, integer *ihi, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    static doublereal aii;
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), xerbla_(char *, integer *, ftnlen), 
	    slarfg_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);


/*  -- LAPACK computational routine (version 3.4.2) -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 184 "sgehd2.f"
    /* Parameter adjustments */
#line 184 "sgehd2.f"
    a_dim1 = *lda;
#line 184 "sgehd2.f"
    a_offset = 1 + a_dim1;
#line 184 "sgehd2.f"
    a -= a_offset;
#line 184 "sgehd2.f"
    --tau;
#line 184 "sgehd2.f"
    --work;
#line 184 "sgehd2.f"

#line 184 "sgehd2.f"
    /* Function Body */
#line 184 "sgehd2.f"
    *info = 0;
#line 185 "sgehd2.f"
    if (*n < 0) {
#line 186 "sgehd2.f"
	*info = -1;
#line 187 "sgehd2.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 188 "sgehd2.f"
	*info = -2;
#line 189 "sgehd2.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 190 "sgehd2.f"
	*info = -3;
#line 191 "sgehd2.f"
    } else if (*lda < max(1,*n)) {
#line 192 "sgehd2.f"
	*info = -5;
#line 193 "sgehd2.f"
    }
#line 194 "sgehd2.f"
    if (*info != 0) {
#line 195 "sgehd2.f"
	i__1 = -(*info);
#line 195 "sgehd2.f"
	xerbla_("SGEHD2", &i__1, (ftnlen)6);
#line 196 "sgehd2.f"
	return 0;
#line 197 "sgehd2.f"
    }

#line 199 "sgehd2.f"
    i__1 = *ihi - 1;
#line 199 "sgehd2.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {

/*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i) */

#line 203 "sgehd2.f"
	i__2 = *ihi - i__;
/* Computing MIN */
#line 203 "sgehd2.f"
	i__3 = i__ + 2;
#line 203 "sgehd2.f"
	slarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ * 
		a_dim1], &c__1, &tau[i__]);
#line 205 "sgehd2.f"
	aii = a[i__ + 1 + i__ * a_dim1];
#line 206 "sgehd2.f"
	a[i__ + 1 + i__ * a_dim1] = 1.;

/*        Apply H(i) to A(1:ihi,i+1:ihi) from the right */

#line 210 "sgehd2.f"
	i__2 = *ihi - i__;
#line 210 "sgehd2.f"
	slarf_("Right", ihi, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
		i__], &a[(i__ + 1) * a_dim1 + 1], lda, &work[1], (ftnlen)5);

/*        Apply H(i) to A(i+1:ihi,i+1:n) from the left */

#line 215 "sgehd2.f"
	i__2 = *ihi - i__;
#line 215 "sgehd2.f"
	i__3 = *n - i__;
#line 215 "sgehd2.f"
	slarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
		i__], &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1], (
		ftnlen)4);

#line 218 "sgehd2.f"
	a[i__ + 1 + i__ * a_dim1] = aii;
#line 219 "sgehd2.f"
/* L10: */
#line 219 "sgehd2.f"
    }

#line 221 "sgehd2.f"
    return 0;

/*     End of SGEHD2 */

} /* sgehd2_ */

