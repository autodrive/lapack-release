#line 1 "dlarf.f"
/* dlarf.f -- translated by f2c (version 20100827).
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

#line 1 "dlarf.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;

/* > \brief \b DLARF applies an elementary reflector to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            INCV, LDC, M, N */
/*       DOUBLE PRECISION   TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARF applies a real elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* >       H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': form  H * C */
/* >          = 'R': form  C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension */
/* >                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
/* >                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
/* >          The vector v in the representation of H. V is not used if */
/* >          TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* >          INCV is INTEGER */
/* >          The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* >          or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension */
/* >                         (N) if SIDE = 'L' */
/* >                      or (M) if SIDE = 'R' */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlarf_(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work, ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static logical applyleft;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer lastc, lastv;
    extern integer iladlc_(integer *, integer *, doublereal *, integer *), 
	    iladlr_(integer *, integer *, doublereal *, integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 161 "dlarf.f"
    /* Parameter adjustments */
#line 161 "dlarf.f"
    --v;
#line 161 "dlarf.f"
    c_dim1 = *ldc;
#line 161 "dlarf.f"
    c_offset = 1 + c_dim1;
#line 161 "dlarf.f"
    c__ -= c_offset;
#line 161 "dlarf.f"
    --work;
#line 161 "dlarf.f"

#line 161 "dlarf.f"
    /* Function Body */
#line 161 "dlarf.f"
    applyleft = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 162 "dlarf.f"
    lastv = 0;
#line 163 "dlarf.f"
    lastc = 0;
#line 164 "dlarf.f"
    if (*tau != 0.) {
/*     Set up variables for scanning V.  LASTV begins pointing to the end */
/*     of V. */
#line 167 "dlarf.f"
	if (applyleft) {
#line 168 "dlarf.f"
	    lastv = *m;
#line 169 "dlarf.f"
	} else {
#line 170 "dlarf.f"
	    lastv = *n;
#line 171 "dlarf.f"
	}
#line 172 "dlarf.f"
	if (*incv > 0) {
#line 173 "dlarf.f"
	    i__ = (lastv - 1) * *incv + 1;
#line 174 "dlarf.f"
	} else {
#line 175 "dlarf.f"
	    i__ = 1;
#line 176 "dlarf.f"
	}
/*     Look for the last non-zero row in V. */
#line 178 "dlarf.f"
	while(lastv > 0 && v[i__] == 0.) {
#line 179 "dlarf.f"
	    --lastv;
#line 180 "dlarf.f"
	    i__ -= *incv;
#line 181 "dlarf.f"
	}
#line 182 "dlarf.f"
	if (applyleft) {
/*     Scan for the last non-zero column in C(1:lastv,:). */
#line 184 "dlarf.f"
	    lastc = iladlc_(&lastv, n, &c__[c_offset], ldc);
#line 185 "dlarf.f"
	} else {
/*     Scan for the last non-zero row in C(:,1:lastv). */
#line 187 "dlarf.f"
	    lastc = iladlr_(m, &lastv, &c__[c_offset], ldc);
#line 188 "dlarf.f"
	}
#line 189 "dlarf.f"
    }
/*     Note that lastc.eq.0 renders the BLAS operations null; no special */
/*     case is needed at this level. */
#line 192 "dlarf.f"
    if (applyleft) {

/*        Form  H * C */

#line 196 "dlarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */

#line 200 "dlarf.f"
	    dgemv_("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &
		    v[1], incv, &c_b5, &work[1], &c__1, (ftnlen)9);

/*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T */

#line 205 "dlarf.f"
	    d__1 = -(*tau);
#line 205 "dlarf.f"
	    dger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[
		    c_offset], ldc);
#line 206 "dlarf.f"
	}
#line 207 "dlarf.f"
    } else {

/*        Form  C * H */

#line 211 "dlarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */

#line 215 "dlarf.f"
	    dgemv_("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc,
		     &v[1], incv, &c_b5, &work[1], &c__1, (ftnlen)12);

/*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */

#line 220 "dlarf.f"
	    d__1 = -(*tau);
#line 220 "dlarf.f"
	    dger_(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[
		    c_offset], ldc);
#line 221 "dlarf.f"
	}
#line 222 "dlarf.f"
    }
#line 223 "dlarf.f"
    return 0;

/*     End of DLARF */

} /* dlarf_ */

