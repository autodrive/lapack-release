#line 1 "zunghr.f"
/* zunghr.f -- translated by f2c (version 20100827).
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

#line 1 "zunghr.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZUNGHR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGHR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunghr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunghr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunghr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNGHR generates a complex unitary matrix Q which is defined as the */
/* > product of IHI-ILO elementary reflectors of order N, as returned by */
/* > ZGEHRD: */
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
/* >          of ZGEHRD. Q is equal to the unit matrix except in the */
/* >          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
/* >          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by ZGEHRD. */
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
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGEHRD. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zunghr_(integer *n, integer *ilo, integer *ihi, 
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
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);


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

#line 166 "zunghr.f"
    /* Parameter adjustments */
#line 166 "zunghr.f"
    a_dim1 = *lda;
#line 166 "zunghr.f"
    a_offset = 1 + a_dim1;
#line 166 "zunghr.f"
    a -= a_offset;
#line 166 "zunghr.f"
    --tau;
#line 166 "zunghr.f"
    --work;
#line 166 "zunghr.f"

#line 166 "zunghr.f"
    /* Function Body */
#line 166 "zunghr.f"
    *info = 0;
#line 167 "zunghr.f"
    nh = *ihi - *ilo;
#line 168 "zunghr.f"
    lquery = *lwork == -1;
#line 169 "zunghr.f"
    if (*n < 0) {
#line 170 "zunghr.f"
	*info = -1;
#line 171 "zunghr.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 172 "zunghr.f"
	*info = -2;
#line 173 "zunghr.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 174 "zunghr.f"
	*info = -3;
#line 175 "zunghr.f"
    } else if (*lda < max(1,*n)) {
#line 176 "zunghr.f"
	*info = -5;
#line 177 "zunghr.f"
    } else if (*lwork < max(1,nh) && ! lquery) {
#line 178 "zunghr.f"
	*info = -8;
#line 179 "zunghr.f"
    }

#line 181 "zunghr.f"
    if (*info == 0) {
#line 182 "zunghr.f"
	nb = ilaenv_(&c__1, "ZUNGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (
		ftnlen)1);
#line 183 "zunghr.f"
	lwkopt = max(1,nh) * nb;
#line 184 "zunghr.f"
	work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 185 "zunghr.f"
    }

#line 187 "zunghr.f"
    if (*info != 0) {
#line 188 "zunghr.f"
	i__1 = -(*info);
#line 188 "zunghr.f"
	xerbla_("ZUNGHR", &i__1, (ftnlen)6);
#line 189 "zunghr.f"
	return 0;
#line 190 "zunghr.f"
    } else if (lquery) {
#line 191 "zunghr.f"
	return 0;
#line 192 "zunghr.f"
    }

/*     Quick return if possible */

#line 196 "zunghr.f"
    if (*n == 0) {
#line 197 "zunghr.f"
	work[1].r = 1., work[1].i = 0.;
#line 198 "zunghr.f"
	return 0;
#line 199 "zunghr.f"
    }

/*     Shift the vectors which define the elementary reflectors one */
/*     column to the right, and set the first ilo and the last n-ihi */
/*     rows and columns to those of the unit matrix */

#line 205 "zunghr.f"
    i__1 = *ilo + 1;
#line 205 "zunghr.f"
    for (j = *ihi; j >= i__1; --j) {
#line 206 "zunghr.f"
	i__2 = j - 1;
#line 206 "zunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 207 "zunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 207 "zunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 208 "zunghr.f"
/* L10: */
#line 208 "zunghr.f"
	}
#line 209 "zunghr.f"
	i__2 = *ihi;
#line 209 "zunghr.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 210 "zunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 210 "zunghr.f"
	    i__4 = i__ + (j - 1) * a_dim1;
#line 210 "zunghr.f"
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
#line 211 "zunghr.f"
/* L20: */
#line 211 "zunghr.f"
	}
#line 212 "zunghr.f"
	i__2 = *n;
#line 212 "zunghr.f"
	for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
#line 213 "zunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 213 "zunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 214 "zunghr.f"
/* L30: */
#line 214 "zunghr.f"
	}
#line 215 "zunghr.f"
/* L40: */
#line 215 "zunghr.f"
    }
#line 216 "zunghr.f"
    i__1 = *ilo;
#line 216 "zunghr.f"
    for (j = 1; j <= i__1; ++j) {
#line 217 "zunghr.f"
	i__2 = *n;
#line 217 "zunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 218 "zunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 218 "zunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 219 "zunghr.f"
/* L50: */
#line 219 "zunghr.f"
	}
#line 220 "zunghr.f"
	i__2 = j + j * a_dim1;
#line 220 "zunghr.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 221 "zunghr.f"
/* L60: */
#line 221 "zunghr.f"
    }
#line 222 "zunghr.f"
    i__1 = *n;
#line 222 "zunghr.f"
    for (j = *ihi + 1; j <= i__1; ++j) {
#line 223 "zunghr.f"
	i__2 = *n;
#line 223 "zunghr.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 224 "zunghr.f"
	    i__3 = i__ + j * a_dim1;
#line 224 "zunghr.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 225 "zunghr.f"
/* L70: */
#line 225 "zunghr.f"
	}
#line 226 "zunghr.f"
	i__2 = j + j * a_dim1;
#line 226 "zunghr.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 227 "zunghr.f"
/* L80: */
#line 227 "zunghr.f"
    }

#line 229 "zunghr.f"
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

#line 233 "zunghr.f"
	zungqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*
		ilo], &work[1], lwork, &iinfo);
#line 235 "zunghr.f"
    }
#line 236 "zunghr.f"
    work[1].r = (doublereal) lwkopt, work[1].i = 0.;
#line 237 "zunghr.f"
    return 0;

/*     End of ZUNGHR */

} /* zunghr_ */

