#line 1 "ctrexc.f"
/* ctrexc.f -- translated by f2c (version 20100827).
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

#line 1 "ctrexc.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTREXC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTREXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrexc.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrexc.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrexc.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPQ */
/*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            Q( LDQ, * ), T( LDT, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTREXC reorders the Schur factorization of a complex matrix */
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
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,N) */
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
/* >          Q is COMPLEX array, dimension (LDQ,N) */
/* >          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
/* >          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
/* >          unitary transformation matrix Z which reorders T. */
/* >          If COMPQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= max(1,N). */
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

/* > \date November 2011 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctrexc_(char *compq, integer *n, doublecomplex *t, 
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
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq;
    extern /* Subroutine */ int clartg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *, doublecomplex *), xerbla_(char *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 162 "ctrexc.f"
    /* Parameter adjustments */
#line 162 "ctrexc.f"
    t_dim1 = *ldt;
#line 162 "ctrexc.f"
    t_offset = 1 + t_dim1;
#line 162 "ctrexc.f"
    t -= t_offset;
#line 162 "ctrexc.f"
    q_dim1 = *ldq;
#line 162 "ctrexc.f"
    q_offset = 1 + q_dim1;
#line 162 "ctrexc.f"
    q -= q_offset;
#line 162 "ctrexc.f"

#line 162 "ctrexc.f"
    /* Function Body */
#line 162 "ctrexc.f"
    *info = 0;
#line 163 "ctrexc.f"
    wantq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 164 "ctrexc.f"
    if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 165 "ctrexc.f"
	*info = -1;
#line 166 "ctrexc.f"
    } else if (*n < 0) {
#line 167 "ctrexc.f"
	*info = -2;
#line 168 "ctrexc.f"
    } else if (*ldt < max(1,*n)) {
#line 169 "ctrexc.f"
	*info = -4;
#line 170 "ctrexc.f"
    } else if (*ldq < 1 || wantq && *ldq < max(1,*n)) {
#line 171 "ctrexc.f"
	*info = -6;
#line 172 "ctrexc.f"
    } else if (*ifst < 1 || *ifst > *n) {
#line 173 "ctrexc.f"
	*info = -7;
#line 174 "ctrexc.f"
    } else if (*ilst < 1 || *ilst > *n) {
#line 175 "ctrexc.f"
	*info = -8;
#line 176 "ctrexc.f"
    }
#line 177 "ctrexc.f"
    if (*info != 0) {
#line 178 "ctrexc.f"
	i__1 = -(*info);
#line 178 "ctrexc.f"
	xerbla_("CTREXC", &i__1, (ftnlen)6);
#line 179 "ctrexc.f"
	return 0;
#line 180 "ctrexc.f"
    }

/*     Quick return if possible */

#line 184 "ctrexc.f"
    if (*n == 1 || *ifst == *ilst) {
#line 184 "ctrexc.f"
	return 0;
#line 184 "ctrexc.f"
    }

#line 187 "ctrexc.f"
    if (*ifst < *ilst) {

/*        Move the IFST-th diagonal element forward down the diagonal. */

#line 191 "ctrexc.f"
	m1 = 0;
#line 192 "ctrexc.f"
	m2 = -1;
#line 193 "ctrexc.f"
	m3 = 1;
#line 194 "ctrexc.f"
    } else {

/*        Move the IFST-th diagonal element backward up the diagonal. */

#line 198 "ctrexc.f"
	m1 = -1;
#line 199 "ctrexc.f"
	m2 = 0;
#line 200 "ctrexc.f"
	m3 = -1;
#line 201 "ctrexc.f"
    }

#line 203 "ctrexc.f"
    i__1 = *ilst + m2;
#line 203 "ctrexc.f"
    i__2 = m3;
#line 203 "ctrexc.f"
    for (k = *ifst + m1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {

/*        Interchange the k-th and (k+1)-th diagonal elements. */

#line 207 "ctrexc.f"
	i__3 = k + k * t_dim1;
#line 207 "ctrexc.f"
	t11.r = t[i__3].r, t11.i = t[i__3].i;
#line 208 "ctrexc.f"
	i__3 = k + 1 + (k + 1) * t_dim1;
#line 208 "ctrexc.f"
	t22.r = t[i__3].r, t22.i = t[i__3].i;

/*        Determine the transformation to perform the interchange. */

#line 212 "ctrexc.f"
	z__1.r = t22.r - t11.r, z__1.i = t22.i - t11.i;
#line 212 "ctrexc.f"
	clartg_(&t[k + (k + 1) * t_dim1], &z__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix T. */

#line 216 "ctrexc.f"
	if (k + 2 <= *n) {
#line 216 "ctrexc.f"
	    i__3 = *n - k - 1;
#line 216 "ctrexc.f"
	    crot_(&i__3, &t[k + (k + 2) * t_dim1], ldt, &t[k + 1 + (k + 2) * 
		    t_dim1], ldt, &cs, &sn);
#line 216 "ctrexc.f"
	}
#line 219 "ctrexc.f"
	i__3 = k - 1;
#line 219 "ctrexc.f"
	d_cnjg(&z__1, &sn);
#line 219 "ctrexc.f"
	crot_(&i__3, &t[k * t_dim1 + 1], &c__1, &t[(k + 1) * t_dim1 + 1], &
		c__1, &cs, &z__1);

#line 221 "ctrexc.f"
	i__3 = k + k * t_dim1;
#line 221 "ctrexc.f"
	t[i__3].r = t22.r, t[i__3].i = t22.i;
#line 222 "ctrexc.f"
	i__3 = k + 1 + (k + 1) * t_dim1;
#line 222 "ctrexc.f"
	t[i__3].r = t11.r, t[i__3].i = t11.i;

#line 224 "ctrexc.f"
	if (wantq) {

/*           Accumulate transformation in the matrix Q. */

#line 228 "ctrexc.f"
	    d_cnjg(&z__1, &sn);
#line 228 "ctrexc.f"
	    crot_(n, &q[k * q_dim1 + 1], &c__1, &q[(k + 1) * q_dim1 + 1], &
		    c__1, &cs, &z__1);
#line 230 "ctrexc.f"
	}

#line 232 "ctrexc.f"
/* L10: */
#line 232 "ctrexc.f"
    }

#line 234 "ctrexc.f"
    return 0;

/*     End of CTREXC */

} /* ctrexc_ */

