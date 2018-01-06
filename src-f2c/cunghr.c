#line 1 "cunghr.f"
/* cunghr.f -- translated by f2c (version 20100827).
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

#line 1 "cunghr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b CUNGHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CUNGHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunghr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunghr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunghr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNGHR generates a complex unitary matrix Q which is defined as the */
/* > product of IHI-ILO elementary reflectors of order N, as returned by */
/* > CGEHRD: */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix Q. N >= 0. */
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
/* >          ILO and IHI must have the same values as in the previous call */
/* >          of CGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by CGEHRD. */
/* >          On exit, the N-by-N unitary matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by CGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= IHI-ILO. */
/* >          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is */
/* >          the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cunghr_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, nb, nh, iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static integer lwkopt;
    static logical lquery;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 166 "cunghr.f"
    /* Parameter adjustments */
#line 166 "cunghr.f"
    a_dim1 = *lda;
#line 166 "cunghr.f"
    a_offset = 1 + a_dim1;
#line 166 "cunghr.f"
    a -= a_offset;
#line 166 "cunghr.f"
    --tau;
#line 166 "cunghr.f"
    --work;
#line 166 "cunghr.f"

#line 166 "cunghr.f"
    /* Function Body */
#line 166 "cunghr.f"
    *info = 0;
#line 167 "cunghr.f"
    nh = *ihi - *ilo;
#line 168 "cunghr.f"
    lquery = *lwork == -1;
#line 169 "cunghr.f"
    if (*n < 0) {
#line 170 "cunghr.f"
	*info = -1;
#line 171 "cunghr.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 172 "cunghr.f"
	*info = -2;
#line 173 "cunghr.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 174 "cunghr.f"
	*info = -3;
#line 175 "cunghr.f"
    } else if (*lda < max(1,*n)) {
#line 176 "cunghr.f"
	*info = -5;
#line 177 "cunghr.f"
    } else if (*lwork < max(1,nh) && ! lquery) {
#line 178 "cunghr.f"
	*info = -8;
#line 179 "cunghr.f"
    }

#line 181 "cunghr.f"
    if (*info == 0) {
#line 182 "cunghr.f"
	nb = ilaenv_(&c__1, "CUNGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 183 "cunghr.f"
	lwkopt = max(1,nh) * nb;
#line 184 "cunghr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 185 "cunghr.f"
    }

#line 187 "cunghr.f"
    if (*info != 0) {
#line 188 "cunghr.f"
	i__1 = -(*info);
#line 188 "cunghr.f"
	xerbla_("CUNGHR", &i__1, (ftnlen)6);
#line 189 "cunghr.f"
	return 0;
#line 190 "cunghr.f"
    } else if (lquery) {
#line 191 "cunghr.f"
	return 0;
#line 192 "cunghr.f"
    }

/*     Quick return if possible */

#line 196 "cunghr.f"
    if (*n == 0) {
#line 197 "cunghr.f"
	work[1].r = 1., work[1].i = 0.;
#line 198 "cunghr.f"
	return 0;
#line 199 "cunghr.f"
    }

/*     Shift the vectors which define the elementary reflectors one */
/*     column to the right, and set the first ilo and the last n-ihi */
/*     rows and columns to those of the unit matrix */

#line 205 "cunghr.f"
    i__1 = *ilo + 1;
#line 205 "cunghr.f"
    for (j = *ihi; j >= i__1; --j) {
#line 206 "cunghr.f"
	i__2 = j - 1;
#line 206 "cunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 207 "cunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 207 "cunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 208 "cunghr.f"
/* L10: */
#line 208 "cunghr.f"
	}
#line 209 "cunghr.f"
	i__2 = *ihi;
#line 209 "cunghr.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 210 "cunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 210 "cunghr.f"
	    i__4 = i__ + (j - 1) * a_dim1;
#line 210 "cunghr.f"
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 211 "cunghr.f"
/* L20: */
#line 211 "cunghr.f"
	}
#line 212 "cunghr.f"
	i__2 = *n;
#line 212 "cunghr.f"
	for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
#line 213 "cunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 213 "cunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 214 "cunghr.f"
/* L30: */
#line 214 "cunghr.f"
	}
#line 215 "cunghr.f"
/* L40: */
#line 215 "cunghr.f"
    }
#line 216 "cunghr.f"
    i__1 = *ilo;
#line 216 "cunghr.f"
    for (j = 1; j <= i__1; ++j) {
#line 217 "cunghr.f"
	i__2 = *n;
#line 217 "cunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 218 "cunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 218 "cunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 219 "cunghr.f"
/* L50: */
#line 219 "cunghr.f"
	}
#line 220 "cunghr.f"
	i__2 = j + j * a_dim1;
#line 220 "cunghr.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 221 "cunghr.f"
/* L60: */
#line 221 "cunghr.f"
    }
#line 222 "cunghr.f"
    i__1 = *n;
#line 222 "cunghr.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 223 "cunghr.f"
	i__2 = *n;
#line 223 "cunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 224 "cunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 224 "cunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 225 "cunghr.f"
/* L70: */
#line 225 "cunghr.f"
	}
#line 226 "cunghr.f"
	i__2 = j + j * a_dim1;
#line 226 "cunghr.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 227 "cunghr.f"
/* L80: */
#line 227 "cunghr.f"
    }

#line 229 "cunghr.f"
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

#line 233 "cunghr.f"
	cungqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*
		ilo], &work[1], lwork, &iinfo);
#line 235 "cunghr.f"
    }
#line 236 "cunghr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 237 "cunghr.f"
    return 0;

/*     End of CUNGHR */

} /* cunghr_ */

