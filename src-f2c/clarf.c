#line 1 "clarf.f"
/* clarf.f -- translated by f2c (version 20100827).
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

#line 1 "clarf.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b CLARF applies an elementary reflector to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            INCV, LDC, M, N */
/*       COMPLEX            TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARF applies a complex elementary reflector H to a complex M-by-N */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* >       H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H (the conjugate transpose of H), supply conjg(tau) instead */
/* > tau. */
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
/* >          V is COMPLEX array, dimension */
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
/* >          TAU is COMPLEX */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
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
/* >          WORK is COMPLEX array, dimension */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clarf_(char *side, integer *m, integer *n, doublecomplex 
	*v, integer *incv, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work, ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;
    doublecomplex z__1;

    /* Local variables */
    static integer i__;
    static logical applyleft;
    extern /* Subroutine */ int cgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer lastc, lastv;
    extern integer ilaclc_(integer *, integer *, doublecomplex *, integer *), 
	    ilaclr_(integer *, integer *, doublecomplex *, integer *);


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

#line 166 "clarf.f"
    /* Parameter adjustments */
#line 166 "clarf.f"
    --v;
#line 166 "clarf.f"
    c_dim1 = *ldc;
#line 166 "clarf.f"
    c_offset = 1 + c_dim1;
#line 166 "clarf.f"
    c__ -= c_offset;
#line 166 "clarf.f"
    --work;
#line 166 "clarf.f"

#line 166 "clarf.f"
    /* Function Body */
#line 166 "clarf.f"
    applyleft = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 167 "clarf.f"
    lastv = 0;
#line 168 "clarf.f"
    lastc = 0;
#line 169 "clarf.f"
    if (tau->r != 0. || tau->i != 0.) {
/*     Set up variables for scanning V.  LASTV begins pointing to the end */
/*     of V. */
#line 172 "clarf.f"
	if (applyleft) {
#line 173 "clarf.f"
	    lastv = *m;
#line 174 "clarf.f"
	} else {
#line 175 "clarf.f"
	    lastv = *n;
#line 176 "clarf.f"
	}
#line 177 "clarf.f"
	if (*incv > 0) {
#line 178 "clarf.f"
	    i__ = (lastv - 1) * *incv + 1;
#line 179 "clarf.f"
	} else {
#line 180 "clarf.f"
	    i__ = 1;
#line 181 "clarf.f"
	}
/*     Look for the last non-zero row in V. */
#line 183 "clarf.f"
	for(;;) { /* while(complicated condition) */
#line 183 "clarf.f"
	    i__1 = i__;
#line 183 "clarf.f"
	    if (!(lastv > 0 && (v[i__1].r == 0. && v[i__1].i == 0.)))
#line 183 "clarf.f"
	    	break;
#line 184 "clarf.f"
	    --lastv;
#line 185 "clarf.f"
	    i__ -= *incv;
#line 186 "clarf.f"
	}
#line 187 "clarf.f"
	if (applyleft) {
/*     Scan for the last non-zero column in C(1:lastv,:). */
#line 189 "clarf.f"
	    lastc = ilaclc_(&lastv, n, &c__[c_offset], ldc);
#line 190 "clarf.f"
	} else {
/*     Scan for the last non-zero row in C(:,1:lastv). */
#line 192 "clarf.f"
	    lastc = ilaclr_(m, &lastv, &c__[c_offset], ldc);
#line 193 "clarf.f"
	}
#line 194 "clarf.f"
    }
/*     Note that lastc.eq.0 renders the BLAS operations null; no special */
/*     case is needed at this level. */
#line 197 "clarf.f"
    if (applyleft) {

/*        Form  H * C */

#line 201 "clarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1) */

#line 205 "clarf.f"
	    cgemv_("Conjugate transpose", &lastv, &lastc, &c_b1, &c__[
		    c_offset], ldc, &v[1], incv, &c_b2, &work[1], &c__1, (
		    ftnlen)19);

/*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H */

#line 210 "clarf.f"
	    z__1.r = -tau->r, z__1.i = -tau->i;
#line 210 "clarf.f"
	    cgerc_(&lastv, &lastc, &z__1, &v[1], incv, &work[1], &c__1, &c__[
		    c_offset], ldc);
#line 211 "clarf.f"
	}
#line 212 "clarf.f"
    } else {

/*        Form  C * H */

#line 216 "clarf.f"
	if (lastv > 0) {

/*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */

#line 220 "clarf.f"
	    cgemv_("No transpose", &lastc, &lastv, &c_b1, &c__[c_offset], ldc,
		     &v[1], incv, &c_b2, &work[1], &c__1, (ftnlen)12);

/*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H */

#line 225 "clarf.f"
	    z__1.r = -tau->r, z__1.i = -tau->i;
#line 225 "clarf.f"
	    cgerc_(&lastc, &lastv, &z__1, &work[1], &c__1, &v[1], incv, &c__[
		    c_offset], ldc);
#line 226 "clarf.f"
	}
#line 227 "clarf.f"
    }
#line 228 "clarf.f"
    return 0;

/*     End of CLARF */

} /* clarf_ */

