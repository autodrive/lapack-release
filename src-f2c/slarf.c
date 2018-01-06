#line 1 "slarf.f"
/* slarf.f -- translated by f2c (version 20100827).
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

#line 1 "slarf.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;

/* > \brief \b SLARF applies an elementary reflector to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarf.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarf.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarf.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            INCV, LDC, M, N */
/*       REAL               TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARF applies a real elementary reflector H to a real m by n matrix */
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
/* >          V is REAL array, dimension */
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
/* >          TAU is REAL */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array, dimension */
/* >                         (N) if SIDE = 'L' */
/* >                      or (M) if SIDE = 'R' */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slarf_(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work, ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static logical applyleft;
    extern /* Subroutine */ int sger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lastc;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer lastv;
    extern integer ilaslc_(integer *, integer *, doublereal *, integer *), 
	    ilaslr_(integer *, integer *, doublereal *, integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 161 "slarf.f"
    /* Parameter adjustments */
#line 161 "slarf.f"
    --v;
#line 161 "slarf.f"
    c_dim1 = *ldc;
#line 161 "slarf.f"
    c_offset = 1 + c_dim1;
#line 161 "slarf.f"
    c__ -= c_offset;
#line 161 "slarf.f"
    --work;
#line 161 "slarf.f"

#line 161 "slarf.f"
    /* Function Body */
#line 161 "slarf.f"
    applyleft = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 162 "slarf.f"
    lastv = 0;
#line 163 "slarf.f"
    lastc = 0;
#line 164 "slarf.f"
    if (*tau != 0.) {
/*     Set up variables for scanning V.  LASTV begins pointing to the end */
/*     of V. */
#line 167 "slarf.f"
	if (applyleft) {
#line 168 "slarf.f"
	    lastv = *m;
#line 169 "slarf.f"
	} else {
#line 170 "slarf.f"
	    lastv = *n;
#line 171 "slarf.f"
	}
#line 172 "slarf.f"
	if (*incv > 0) {
#line 173 "slarf.f"
	    i__ = (lastv - 1) * *incv + 1;
#line 174 "slarf.f"
	} else {
#line 175 "slarf.f"
	    i__ = 1;
#line 176 "slarf.f"
	}
/*     Look for the last non-zero row in V. */
#line 178 "slarf.f"
	while(lastv > 0 && v[i__] == 0.) {
#line 179 "slarf.f"
	    --lastv;
#line 180 "slarf.f"
	    i__ -= *incv;
#line 181 "slarf.f"
	}
#line 182 "slarf.f"
	if (applyleft) {
/*     Scan for the last non-zero column in C(1:lastv,:). */
#line 184 "slarf.f"
	    lastc = ilaslc_(&lastv, n, &c__[c_offset], ldc);
#line 185 "slarf.f"
	} else {
/*     Scan for the last non-zero row in C(:,1:lastv). */
#line 187 "slarf.f"
	    lastc = ilaslr_(m, &lastv, &c__[c_offset], ldc);
#line 188 "slarf.f"
	}
#line 189 "slarf.f"
    }
/*     Note that lastc.eq.0 renders the BLAS operations null; no special */
/*     case is needed at this level. */
#line 192 "slarf.f"
    if (applyleft) {

/*        Form  H * C */

#line 196 "slarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1) */

#line 200 "slarf.f"
	    sgemv_("Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &
		    v[1], incv, &c_b5, &work[1], &c__1, (ftnlen)9);

/*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T */

#line 205 "slarf.f"
	    d__1 = -(*tau);
#line 205 "slarf.f"
	    sger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[
		    c_offset], ldc);
#line 206 "slarf.f"
	}
#line 207 "slarf.f"
    } else {

/*        Form  C * H */

#line 211 "slarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */

#line 215 "slarf.f"
	    sgemv_("No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc,
		     &v[1], incv, &c_b5, &work[1], &c__1, (ftnlen)12);

/*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T */

#line 220 "slarf.f"
	    d__1 = -(*tau);
#line 220 "slarf.f"
	    sger_(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[
		    c_offset], ldc);
#line 221 "slarf.f"
	}
#line 222 "slarf.f"
    }
#line 223 "slarf.f"
    return 0;

/*     End of SLARF */

} /* slarf_ */

