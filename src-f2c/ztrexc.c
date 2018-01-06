#line 1 "ztrexc.f"
/* ztrexc.f -- translated by f2c (version 20100827).
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

#line 1 "ztrexc.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTREXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTREXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ */
/*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         Q( LDQ, * ), T( LDT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTREXC reorders the Schur factorization of a complex matrix */
/* > A = Q*T*Q**H, so that the diagonal element of T with row index IFST */
/* > is moved to row ILST. */
/* > */
/* > The Schur form T is reordered by a unitary similarity transformation */
/* > Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by */
/* > postmultplying it with Z. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPQ */
/* > \verbatim */
/* >          COMPQ is CHARACTER*1 */
/* >          = 'V':  update the matrix Q of Schur vectors; */
/* >          = 'N':  do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix T. N >= 0. */
/* >          If N == 0 arguments ILST and IFST may be any value. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,N) */
/* >          On entry, the upper triangular matrix T. */
/* >          On exit, the reordered upper triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
/* >          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
/* >          unitary transformation matrix Z which reorders T. */
/* >          If COMPQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= 1, and if */
/* >          COMPQ = 'V', LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IFST */
/* > \verbatim */
/* >          IFST is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] ILST */
/* > \verbatim */
/* >          ILST is INTEGER */
/* > */
/* >          Specify the reordering of the diagonal elements of T: */
/* >          The element with row index IFST is moved to row ILST by a */
/* >          sequence of transpositions between adjacent elements. */
/* >          1 <= IFST <= N; 1 <= ILST <= N. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztrexc_(char *compq, integer *n, doublecomplex *t, 
	integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *
	ilst, integer *info, ftnlen compq_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, t_dim1, t_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, m1, m2, m3;
    static doublereal cs;
    static doublecomplex t11, t22, sn, temp;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlartg_(
	    doublecomplex *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters. */

#line 164 "ztrexc.f"
    /* Parameter adjustments */
#line 164 "ztrexc.f"
    t_dim1 = *ldt;
#line 164 "ztrexc.f"
    t_offset = 1 + t_dim1;
#line 164 "ztrexc.f"
    t -= t_offset;
#line 164 "ztrexc.f"
    q_dim1 = *ldq;
#line 164 "ztrexc.f"
    q_offset = 1 + q_dim1;
#line 164 "ztrexc.f"
    q -= q_offset;
#line 164 "ztrexc.f"

#line 164 "ztrexc.f"
    /* Function Body */
#line 164 "ztrexc.f"
    *info = 0;
#line 165 "ztrexc.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 166 "ztrexc.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 167 "ztrexc.f"
	*info = -1;
#line 168 "ztrexc.f"
    } else if (*n < 0) {
#line 169 "ztrexc.f"
	*info = -2;
#line 170 "ztrexc.f"
    } else if (*ldt < max(1,*n)) {
#line 171 "ztrexc.f"
	*info = -4;
#line 172 "ztrexc.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*n)) {
#line 173 "ztrexc.f"
	*info = -6;
#line 174 "ztrexc.f"
    } else if ((*ifst < 1 || *ifst > *n) && *n > 0) {
#line 175 "ztrexc.f"
	*info = -7;
#line 176 "ztrexc.f"
    } else if ((*ilst < 1 || *ilst > *n) && *n > 0) {
#line 177 "ztrexc.f"
	*info = -8;
#line 178 "ztrexc.f"
    }
#line 179 "ztrexc.f"
    if (*info != 0) {
#line 180 "ztrexc.f"
	i__1 = -(*info);
#line 180 "ztrexc.f"
	xerbla_("ZTREXC", &i__1, (ftnlen)6);
#line 181 "ztrexc.f"
	return 0;
#line 182 "ztrexc.f"
    }

/*     Quick return if possible */

#line 186 "ztrexc.f"
    if (*n <= 1 || *ifst == *ilst) {
#line 186 "ztrexc.f"
	return 0;
#line 186 "ztrexc.f"
    }

#line 189 "ztrexc.f"
    if (*ifst < *ilst) {

/*        Move the IFST-th diagonal element forward down the diagonal. */

#line 193 "ztrexc.f"
	m1 = 0;
#line 194 "ztrexc.f"
	m2 = -1;
#line 195 "ztrexc.f"
	m3 = 1;
#line 196 "ztrexc.f"
    } else {

/*        Move the IFST-th diagonal element backward up the diagonal. */

#line 200 "ztrexc.f"
	m1 = -1;
#line 201 "ztrexc.f"
	m2 = 0;
#line 202 "ztrexc.f"
	m3 = -1;
#line 203 "ztrexc.f"
    }

#line 205 "ztrexc.f"
    i__1 = *ilst + m2;
#line 205 "ztrexc.f"
    i__2 = m3;
#line 205 "ztrexc.f"
    for (k = *ifst + m1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {

/*        Interchange the k-th and (k+1)-th diagonal elements. */

#line 209 "ztrexc.f"
	i__3 = k + k * t_dim1;
#line 209 "ztrexc.f"
	t11.r = t[i__3].r, t11.i = t[i__3].i;
#line 210 "ztrexc.f"
	i__3 = k + 1 + (k + 1) * t_dim1;
#line 210 "ztrexc.f"
	t22.r = t[i__3].r, t22.i = t[i__3].i;

/*        Determine the transformation to perform the interchange. */

#line 214 "ztrexc.f"
	z__1.r = t22.r - t11.r, z__1.i = t22.i - t11.i;
#line 214 "ztrexc.f"
	zlartg_(&t[k + (k + 1) * t_dim1], &z__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix T. */

#line 218 "ztrexc.f"
	if (k + 2 <= *n) {
#line 218 "ztrexc.f"
	    i__3 = *n - k - 1;
#line 218 "ztrexc.f"
	    zrot_(&i__3, &t[k + (k + 2) * t_dim1], ldt, &t[k + 1 + (k + 2) * 
		    t_dim1], ldt, &cs, &sn);
#line 218 "ztrexc.f"
	}
#line 221 "ztrexc.f"
	i__3 = k - 1;
#line 221 "ztrexc.f"
	d_cnjg(&z__1, &sn);
#line 221 "ztrexc.f"
	zrot_(&i__3, &t[k * t_dim1 + 1], &c__1, &t[(k + 1) * t_dim1 + 1], &
		c__1, &cs, &z__1);

#line 224 "ztrexc.f"
	i__3 = k + k * t_dim1;
#line 224 "ztrexc.f"
	t[i__3].r = t22.r, t[i__3].i = t22.i;
#line 225 "ztrexc.f"
	i__3 = k + 1 + (k + 1) * t_dim1;
#line 225 "ztrexc.f"
	t[i__3].r = t11.r, t[i__3].i = t11.i;

#line 227 "ztrexc.f"
	if (wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 231 "ztrexc.f"
	    d_cnjg(&z__1, &sn);
#line 231 "ztrexc.f"
	    zrot_(n, &q[k * q_dim1 + 1], &c__1, &q[(k + 1) * q_dim1 + 1], &
		    c__1, &cs, &z__1);
#line 233 "ztrexc.f"
	}

#line 235 "ztrexc.f"
/* L10: */
#line 235 "ztrexc.f"
    }

#line 237 "ztrexc.f"
    return 0;

/*     End of ZTREXC */

} /* ztrexc_ */

