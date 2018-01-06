#line 1 "zlags2.f"
/* zlags2.f -- translated by f2c (version 20100827).
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

#line 1 "zlags2.f"
/* > \brief \b ZLAGS2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlags2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlags2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlags2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, */
/*                          SNV, CSQ, SNQ ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            UPPER */
/*       DOUBLE PRECISION   A1, A3, B1, B3, CSQ, CSU, CSV */
/*       COMPLEX*16         A2, B2, SNQ, SNU, SNV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAGS2 computes 2-by-2 unitary matrices U, V and Q, such */
/* > that if ( UPPER ) then */
/* > */
/* >           U**H *A*Q = U**H *( A1 A2 )*Q = ( x  0  ) */
/* >                             ( 0  A3 )     ( x  x  ) */
/* > and */
/* >           V**H*B*Q = V**H *( B1 B2 )*Q = ( x  0  ) */
/* >                            ( 0  B3 )     ( x  x  ) */
/* > */
/* > or if ( .NOT.UPPER ) then */
/* > */
/* >           U**H *A*Q = U**H *( A1 0  )*Q = ( x  x  ) */
/* >                             ( A2 A3 )     ( 0  x  ) */
/* > and */
/* >           V**H *B*Q = V**H *( B1 0  )*Q = ( x  x  ) */
/* >                             ( B2 B3 )     ( 0  x  ) */
/* > where */
/* > */
/* >   U = (   CSU    SNU ), V = (  CSV    SNV ), */
/* >       ( -SNU**H  CSU )      ( -SNV**H CSV ) */
/* > */
/* >   Q = (   CSQ    SNQ ) */
/* >       ( -SNQ**H  CSQ ) */
/* > */
/* > The rows of the transformed A and B are parallel. Moreover, if the */
/* > input 2-by-2 matrix A is not zero, then the transformed (1,1) entry */
/* > of A is not zero. If the input matrices A and B are both not zero, */
/* > then the transformed (2,2) element of B is not zero, except when the */
/* > first rows of input A and B are parallel and the second rows are */
/* > zero. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPPER */
/* > \verbatim */
/* >          UPPER is LOGICAL */
/* >          = .TRUE.: the input matrices A and B are upper triangular. */
/* >          = .FALSE.: the input matrices A and B are lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] A1 */
/* > \verbatim */
/* >          A1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] A2 */
/* > \verbatim */
/* >          A2 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] A3 */
/* > \verbatim */
/* >          A3 is DOUBLE PRECISION */
/* >          On entry, A1, A2 and A3 are elements of the input 2-by-2 */
/* >          upper (lower) triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] B1 */
/* > \verbatim */
/* >          B1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] B2 */
/* > \verbatim */
/* >          B2 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] B3 */
/* > \verbatim */
/* >          B3 is DOUBLE PRECISION */
/* >          On entry, B1, B2 and B3 are elements of the input 2-by-2 */
/* >          upper (lower) triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] CSU */
/* > \verbatim */
/* >          CSU is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SNU */
/* > \verbatim */
/* >          SNU is COMPLEX*16 */
/* >          The desired unitary matrix U. */
/* > \endverbatim */
/* > */
/* > \param[out] CSV */
/* > \verbatim */
/* >          CSV is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SNV */
/* > \verbatim */
/* >          SNV is COMPLEX*16 */
/* >          The desired unitary matrix V. */
/* > \endverbatim */
/* > */
/* > \param[out] CSQ */
/* > \verbatim */
/* >          CSQ is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SNQ */
/* > \verbatim */
/* >          SNQ is COMPLEX*16 */
/* >          The desired unitary matrix Q. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlags2_(logical *upper, doublereal *a1, doublecomplex *
	a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3,
	 doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *
	snv, doublereal *csq, doublecomplex *snq)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double z_abs(doublecomplex *), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal a;
    static doublecomplex b, c__;
    static doublereal d__;
    static doublecomplex r__, d1;
    static doublereal s1, s2, fb, fc;
    static doublecomplex ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22;
    static doublereal csl, csr, snl, snr, aua11, aua12, aua21, aua22, avb12, 
	    avb11, avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), zlartg_(doublecomplex *
	    , doublecomplex *, doublereal *, doublecomplex *, doublecomplex *)
	    ;


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 199 "zlags2.f"
    if (*upper) {

/*        Input matrices A and B are upper triangular matrices */

/*        Form matrix C = A*adj(B) = ( a b ) */
/*                                   ( 0 d ) */

#line 206 "zlags2.f"
	a = *a1 * *b3;
#line 207 "zlags2.f"
	d__ = *a3 * *b1;
#line 208 "zlags2.f"
	z__2.r = *b1 * a2->r, z__2.i = *b1 * a2->i;
#line 208 "zlags2.f"
	z__3.r = *a1 * b2->r, z__3.i = *a1 * b2->i;
#line 208 "zlags2.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 208 "zlags2.f"
	b.r = z__1.r, b.i = z__1.i;
#line 209 "zlags2.f"
	fb = z_abs(&b);

/*        Transform complex 2-by-2 matrix C to real matrix by unitary */
/*        diagonal matrix diag(1,D1). */

#line 214 "zlags2.f"
	d1.r = 1., d1.i = 0.;
#line 215 "zlags2.f"
	if (fb != 0.) {
#line 215 "zlags2.f"
	    z__1.r = b.r / fb, z__1.i = b.i / fb;
#line 215 "zlags2.f"
	    d1.r = z__1.r, d1.i = z__1.i;
#line 215 "zlags2.f"
	}

/*        The SVD of real 2 by 2 triangular C */

/*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T ) */

#line 223 "zlags2.f"
	dlasv2_(&a, &fb, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

#line 225 "zlags2.f"
	if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {

/*           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B, */
/*           and (1,2) element of |U|**H *|A| and |V|**H *|B|. */

#line 231 "zlags2.f"
	    ua11r = csl * *a1;
#line 232 "zlags2.f"
	    z__2.r = csl * a2->r, z__2.i = csl * a2->i;
#line 232 "zlags2.f"
	    z__4.r = snl * d1.r, z__4.i = snl * d1.i;
#line 232 "zlags2.f"
	    z__3.r = *a3 * z__4.r, z__3.i = *a3 * z__4.i;
#line 232 "zlags2.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 232 "zlags2.f"
	    ua12.r = z__1.r, ua12.i = z__1.i;

#line 234 "zlags2.f"
	    vb11r = csr * *b1;
#line 235 "zlags2.f"
	    z__2.r = csr * b2->r, z__2.i = csr * b2->i;
#line 235 "zlags2.f"
	    z__4.r = snr * d1.r, z__4.i = snr * d1.i;
#line 235 "zlags2.f"
	    z__3.r = *b3 * z__4.r, z__3.i = *b3 * z__4.i;
#line 235 "zlags2.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 235 "zlags2.f"
	    vb12.r = z__1.r, vb12.i = z__1.i;

#line 237 "zlags2.f"
	    aua12 = abs(csl) * ((d__1 = a2->r, abs(d__1)) + (d__2 = d_imag(a2)
		    , abs(d__2))) + abs(snl) * abs(*a3);
#line 238 "zlags2.f"
	    avb12 = abs(csr) * ((d__1 = b2->r, abs(d__1)) + (d__2 = d_imag(b2)
		    , abs(d__2))) + abs(snr) * abs(*b3);

/*           zero (1,2) elements of U**H *A and V**H *B */

#line 242 "zlags2.f"
	    if (abs(ua11r) + ((d__1 = ua12.r, abs(d__1)) + (d__2 = d_imag(&
		    ua12), abs(d__2))) == 0.) {
#line 243 "zlags2.f"
		z__2.r = vb11r, z__2.i = 0.;
#line 243 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 243 "zlags2.f"
		d_cnjg(&z__3, &vb12);
#line 243 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 245 "zlags2.f"
	    } else if (abs(vb11r) + ((d__1 = vb12.r, abs(d__1)) + (d__2 = 
		    d_imag(&vb12), abs(d__2))) == 0.) {
#line 246 "zlags2.f"
		z__2.r = ua11r, z__2.i = 0.;
#line 246 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 246 "zlags2.f"
		d_cnjg(&z__3, &ua12);
#line 246 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 248 "zlags2.f"
	    } else if (aua12 / (abs(ua11r) + ((d__1 = ua12.r, abs(d__1)) + (
		    d__2 = d_imag(&ua12), abs(d__2)))) <= avb12 / (abs(vb11r) 
		    + ((d__3 = vb12.r, abs(d__3)) + (d__4 = d_imag(&vb12), 
		    abs(d__4))))) {
#line 250 "zlags2.f"
		z__2.r = ua11r, z__2.i = 0.;
#line 250 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 250 "zlags2.f"
		d_cnjg(&z__3, &ua12);
#line 250 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 252 "zlags2.f"
	    } else {
#line 253 "zlags2.f"
		z__2.r = vb11r, z__2.i = 0.;
#line 253 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 253 "zlags2.f"
		d_cnjg(&z__3, &vb12);
#line 253 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 255 "zlags2.f"
	    }

#line 257 "zlags2.f"
	    *csu = csl;
#line 258 "zlags2.f"
	    z__2.r = -d1.r, z__2.i = -d1.i;
#line 258 "zlags2.f"
	    z__1.r = snl * z__2.r, z__1.i = snl * z__2.i;
#line 258 "zlags2.f"
	    snu->r = z__1.r, snu->i = z__1.i;
#line 259 "zlags2.f"
	    *csv = csr;
#line 260 "zlags2.f"
	    z__2.r = -d1.r, z__2.i = -d1.i;
#line 260 "zlags2.f"
	    z__1.r = snr * z__2.r, z__1.i = snr * z__2.i;
#line 260 "zlags2.f"
	    snv->r = z__1.r, snv->i = z__1.i;

#line 262 "zlags2.f"
	} else {

/*           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B, */
/*           and (2,2) element of |U|**H *|A| and |V|**H *|B|. */

#line 267 "zlags2.f"
	    d_cnjg(&z__4, &d1);
#line 267 "zlags2.f"
	    z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 267 "zlags2.f"
	    z__2.r = snl * z__3.r, z__2.i = snl * z__3.i;
#line 267 "zlags2.f"
	    z__1.r = *a1 * z__2.r, z__1.i = *a1 * z__2.i;
#line 267 "zlags2.f"
	    ua21.r = z__1.r, ua21.i = z__1.i;
#line 268 "zlags2.f"
	    d_cnjg(&z__5, &d1);
#line 268 "zlags2.f"
	    z__4.r = -z__5.r, z__4.i = -z__5.i;
#line 268 "zlags2.f"
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
#line 268 "zlags2.f"
	    z__2.r = z__3.r * a2->r - z__3.i * a2->i, z__2.i = z__3.r * a2->i 
		    + z__3.i * a2->r;
#line 268 "zlags2.f"
	    d__1 = csl * *a3;
#line 268 "zlags2.f"
	    z__1.r = z__2.r + d__1, z__1.i = z__2.i;
#line 268 "zlags2.f"
	    ua22.r = z__1.r, ua22.i = z__1.i;

#line 270 "zlags2.f"
	    d_cnjg(&z__4, &d1);
#line 270 "zlags2.f"
	    z__3.r = -z__4.r, z__3.i = -z__4.i;
#line 270 "zlags2.f"
	    z__2.r = snr * z__3.r, z__2.i = snr * z__3.i;
#line 270 "zlags2.f"
	    z__1.r = *b1 * z__2.r, z__1.i = *b1 * z__2.i;
#line 270 "zlags2.f"
	    vb21.r = z__1.r, vb21.i = z__1.i;
#line 271 "zlags2.f"
	    d_cnjg(&z__5, &d1);
#line 271 "zlags2.f"
	    z__4.r = -z__5.r, z__4.i = -z__5.i;
#line 271 "zlags2.f"
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
#line 271 "zlags2.f"
	    z__2.r = z__3.r * b2->r - z__3.i * b2->i, z__2.i = z__3.r * b2->i 
		    + z__3.i * b2->r;
#line 271 "zlags2.f"
	    d__1 = csr * *b3;
#line 271 "zlags2.f"
	    z__1.r = z__2.r + d__1, z__1.i = z__2.i;
#line 271 "zlags2.f"
	    vb22.r = z__1.r, vb22.i = z__1.i;

#line 273 "zlags2.f"
	    aua22 = abs(snl) * ((d__1 = a2->r, abs(d__1)) + (d__2 = d_imag(a2)
		    , abs(d__2))) + abs(csl) * abs(*a3);
#line 274 "zlags2.f"
	    avb22 = abs(snr) * ((d__1 = b2->r, abs(d__1)) + (d__2 = d_imag(b2)
		    , abs(d__2))) + abs(csr) * abs(*b3);

/*           zero (2,2) elements of U**H *A and V**H *B, and then swap. */

#line 278 "zlags2.f"
	    if ((d__1 = ua21.r, abs(d__1)) + (d__2 = d_imag(&ua21), abs(d__2))
		     + ((d__3 = ua22.r, abs(d__3)) + (d__4 = d_imag(&ua22), 
		    abs(d__4))) == 0.) {
#line 279 "zlags2.f"
		d_cnjg(&z__2, &vb21);
#line 279 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 279 "zlags2.f"
		d_cnjg(&z__3, &vb22);
#line 279 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 281 "zlags2.f"
	    } else if ((d__1 = vb21.r, abs(d__1)) + (d__2 = d_imag(&vb21), 
		    abs(d__2)) + z_abs(&vb22) == 0.) {
#line 282 "zlags2.f"
		d_cnjg(&z__2, &ua21);
#line 282 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 282 "zlags2.f"
		d_cnjg(&z__3, &ua22);
#line 282 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 284 "zlags2.f"
	    } else if (aua22 / ((d__1 = ua21.r, abs(d__1)) + (d__2 = d_imag(&
		    ua21), abs(d__2)) + ((d__3 = ua22.r, abs(d__3)) + (d__4 = 
		    d_imag(&ua22), abs(d__4)))) <= avb22 / ((d__5 = vb21.r, 
		    abs(d__5)) + (d__6 = d_imag(&vb21), abs(d__6)) + ((d__7 = 
		    vb22.r, abs(d__7)) + (d__8 = d_imag(&vb22), abs(d__8))))) 
		    {
#line 286 "zlags2.f"
		d_cnjg(&z__2, &ua21);
#line 286 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 286 "zlags2.f"
		d_cnjg(&z__3, &ua22);
#line 286 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 288 "zlags2.f"
	    } else {
#line 289 "zlags2.f"
		d_cnjg(&z__2, &vb21);
#line 289 "zlags2.f"
		z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 289 "zlags2.f"
		d_cnjg(&z__3, &vb22);
#line 289 "zlags2.f"
		zlartg_(&z__1, &z__3, csq, snq, &r__);
#line 291 "zlags2.f"
	    }

#line 293 "zlags2.f"
	    *csu = snl;
#line 294 "zlags2.f"
	    z__1.r = csl * d1.r, z__1.i = csl * d1.i;
#line 294 "zlags2.f"
	    snu->r = z__1.r, snu->i = z__1.i;
#line 295 "zlags2.f"
	    *csv = snr;
#line 296 "zlags2.f"
	    z__1.r = csr * d1.r, z__1.i = csr * d1.i;
#line 296 "zlags2.f"
	    snv->r = z__1.r, snv->i = z__1.i;

#line 298 "zlags2.f"
	}

#line 300 "zlags2.f"
    } else {

/*        Input matrices A and B are lower triangular matrices */

/*        Form matrix C = A*adj(B) = ( a 0 ) */
/*                                   ( c d ) */

#line 307 "zlags2.f"
	a = *a1 * *b3;
#line 308 "zlags2.f"
	d__ = *a3 * *b1;
#line 309 "zlags2.f"
	z__2.r = *b3 * a2->r, z__2.i = *b3 * a2->i;
#line 309 "zlags2.f"
	z__3.r = *a3 * b2->r, z__3.i = *a3 * b2->i;
#line 309 "zlags2.f"
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
#line 309 "zlags2.f"
	c__.r = z__1.r, c__.i = z__1.i;
#line 310 "zlags2.f"
	fc = z_abs(&c__);

/*        Transform complex 2-by-2 matrix C to real matrix by unitary */
/*        diagonal matrix diag(d1,1). */

#line 315 "zlags2.f"
	d1.r = 1., d1.i = 0.;
#line 316 "zlags2.f"
	if (fc != 0.) {
#line 316 "zlags2.f"
	    z__1.r = c__.r / fc, z__1.i = c__.i / fc;
#line 316 "zlags2.f"
	    d1.r = z__1.r, d1.i = z__1.i;
#line 316 "zlags2.f"
	}

/*        The SVD of real 2 by 2 triangular C */

/*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T ) */

#line 324 "zlags2.f"
	dlasv2_(&a, &fc, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

#line 326 "zlags2.f"
	if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {

/*           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B, */
/*           and (2,1) element of |U|**H *|A| and |V|**H *|B|. */

#line 332 "zlags2.f"
	    z__4.r = -d1.r, z__4.i = -d1.i;
#line 332 "zlags2.f"
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
#line 332 "zlags2.f"
	    z__2.r = *a1 * z__3.r, z__2.i = *a1 * z__3.i;
#line 332 "zlags2.f"
	    z__5.r = csr * a2->r, z__5.i = csr * a2->i;
#line 332 "zlags2.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 332 "zlags2.f"
	    ua21.r = z__1.r, ua21.i = z__1.i;
#line 333 "zlags2.f"
	    ua22r = csr * *a3;

#line 335 "zlags2.f"
	    z__4.r = -d1.r, z__4.i = -d1.i;
#line 335 "zlags2.f"
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
#line 335 "zlags2.f"
	    z__2.r = *b1 * z__3.r, z__2.i = *b1 * z__3.i;
#line 335 "zlags2.f"
	    z__5.r = csl * b2->r, z__5.i = csl * b2->i;
#line 335 "zlags2.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 335 "zlags2.f"
	    vb21.r = z__1.r, vb21.i = z__1.i;
#line 336 "zlags2.f"
	    vb22r = csl * *b3;

#line 338 "zlags2.f"
	    aua21 = abs(snr) * abs(*a1) + abs(csr) * ((d__1 = a2->r, abs(d__1)
		    ) + (d__2 = d_imag(a2), abs(d__2)));
#line 339 "zlags2.f"
	    avb21 = abs(snl) * abs(*b1) + abs(csl) * ((d__1 = b2->r, abs(d__1)
		    ) + (d__2 = d_imag(b2), abs(d__2)));

/*           zero (2,1) elements of U**H *A and V**H *B. */

#line 343 "zlags2.f"
	    if ((d__1 = ua21.r, abs(d__1)) + (d__2 = d_imag(&ua21), abs(d__2))
		     + abs(ua22r) == 0.) {
#line 344 "zlags2.f"
		z__1.r = vb22r, z__1.i = 0.;
#line 344 "zlags2.f"
		zlartg_(&z__1, &vb21, csq, snq, &r__);
#line 345 "zlags2.f"
	    } else if ((d__1 = vb21.r, abs(d__1)) + (d__2 = d_imag(&vb21), 
		    abs(d__2)) + abs(vb22r) == 0.) {
#line 346 "zlags2.f"
		z__1.r = ua22r, z__1.i = 0.;
#line 346 "zlags2.f"
		zlartg_(&z__1, &ua21, csq, snq, &r__);
#line 347 "zlags2.f"
	    } else if (aua21 / ((d__1 = ua21.r, abs(d__1)) + (d__2 = d_imag(&
		    ua21), abs(d__2)) + abs(ua22r)) <= avb21 / ((d__3 = 
		    vb21.r, abs(d__3)) + (d__4 = d_imag(&vb21), abs(d__4)) + 
		    abs(vb22r))) {
#line 349 "zlags2.f"
		z__1.r = ua22r, z__1.i = 0.;
#line 349 "zlags2.f"
		zlartg_(&z__1, &ua21, csq, snq, &r__);
#line 350 "zlags2.f"
	    } else {
#line 351 "zlags2.f"
		z__1.r = vb22r, z__1.i = 0.;
#line 351 "zlags2.f"
		zlartg_(&z__1, &vb21, csq, snq, &r__);
#line 352 "zlags2.f"
	    }

#line 354 "zlags2.f"
	    *csu = csr;
#line 355 "zlags2.f"
	    d_cnjg(&z__3, &d1);
#line 355 "zlags2.f"
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 355 "zlags2.f"
	    z__1.r = snr * z__2.r, z__1.i = snr * z__2.i;
#line 355 "zlags2.f"
	    snu->r = z__1.r, snu->i = z__1.i;
#line 356 "zlags2.f"
	    *csv = csl;
#line 357 "zlags2.f"
	    d_cnjg(&z__3, &d1);
#line 357 "zlags2.f"
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
#line 357 "zlags2.f"
	    z__1.r = snl * z__2.r, z__1.i = snl * z__2.i;
#line 357 "zlags2.f"
	    snv->r = z__1.r, snv->i = z__1.i;

#line 359 "zlags2.f"
	} else {

/*           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B, */
/*           and (1,1) element of |U|**H *|A| and |V|**H *|B|. */

#line 364 "zlags2.f"
	    d__1 = csr * *a1;
#line 364 "zlags2.f"
	    d_cnjg(&z__4, &d1);
#line 364 "zlags2.f"
	    z__3.r = snr * z__4.r, z__3.i = snr * z__4.i;
#line 364 "zlags2.f"
	    z__2.r = z__3.r * a2->r - z__3.i * a2->i, z__2.i = z__3.r * a2->i 
		    + z__3.i * a2->r;
#line 364 "zlags2.f"
	    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
#line 364 "zlags2.f"
	    ua11.r = z__1.r, ua11.i = z__1.i;
#line 365 "zlags2.f"
	    d_cnjg(&z__3, &d1);
#line 365 "zlags2.f"
	    z__2.r = snr * z__3.r, z__2.i = snr * z__3.i;
#line 365 "zlags2.f"
	    z__1.r = *a3 * z__2.r, z__1.i = *a3 * z__2.i;
#line 365 "zlags2.f"
	    ua12.r = z__1.r, ua12.i = z__1.i;

#line 367 "zlags2.f"
	    d__1 = csl * *b1;
#line 367 "zlags2.f"
	    d_cnjg(&z__4, &d1);
#line 367 "zlags2.f"
	    z__3.r = snl * z__4.r, z__3.i = snl * z__4.i;
#line 367 "zlags2.f"
	    z__2.r = z__3.r * b2->r - z__3.i * b2->i, z__2.i = z__3.r * b2->i 
		    + z__3.i * b2->r;
#line 367 "zlags2.f"
	    z__1.r = d__1 + z__2.r, z__1.i = z__2.i;
#line 367 "zlags2.f"
	    vb11.r = z__1.r, vb11.i = z__1.i;
#line 368 "zlags2.f"
	    d_cnjg(&z__3, &d1);
#line 368 "zlags2.f"
	    z__2.r = snl * z__3.r, z__2.i = snl * z__3.i;
#line 368 "zlags2.f"
	    z__1.r = *b3 * z__2.r, z__1.i = *b3 * z__2.i;
#line 368 "zlags2.f"
	    vb12.r = z__1.r, vb12.i = z__1.i;

#line 370 "zlags2.f"
	    aua11 = abs(csr) * abs(*a1) + abs(snr) * ((d__1 = a2->r, abs(d__1)
		    ) + (d__2 = d_imag(a2), abs(d__2)));
#line 371 "zlags2.f"
	    avb11 = abs(csl) * abs(*b1) + abs(snl) * ((d__1 = b2->r, abs(d__1)
		    ) + (d__2 = d_imag(b2), abs(d__2)));

/*           zero (1,1) elements of U**H *A and V**H *B, and then swap. */

#line 375 "zlags2.f"
	    if ((d__1 = ua11.r, abs(d__1)) + (d__2 = d_imag(&ua11), abs(d__2))
		     + ((d__3 = ua12.r, abs(d__3)) + (d__4 = d_imag(&ua12), 
		    abs(d__4))) == 0.) {
#line 376 "zlags2.f"
		zlartg_(&vb12, &vb11, csq, snq, &r__);
#line 377 "zlags2.f"
	    } else if ((d__1 = vb11.r, abs(d__1)) + (d__2 = d_imag(&vb11), 
		    abs(d__2)) + ((d__3 = vb12.r, abs(d__3)) + (d__4 = d_imag(
		    &vb12), abs(d__4))) == 0.) {
#line 378 "zlags2.f"
		zlartg_(&ua12, &ua11, csq, snq, &r__);
#line 379 "zlags2.f"
	    } else if (aua11 / ((d__1 = ua11.r, abs(d__1)) + (d__2 = d_imag(&
		    ua11), abs(d__2)) + ((d__3 = ua12.r, abs(d__3)) + (d__4 = 
		    d_imag(&ua12), abs(d__4)))) <= avb11 / ((d__5 = vb11.r, 
		    abs(d__5)) + (d__6 = d_imag(&vb11), abs(d__6)) + ((d__7 = 
		    vb12.r, abs(d__7)) + (d__8 = d_imag(&vb12), abs(d__8))))) 
		    {
#line 381 "zlags2.f"
		zlartg_(&ua12, &ua11, csq, snq, &r__);
#line 382 "zlags2.f"
	    } else {
#line 383 "zlags2.f"
		zlartg_(&vb12, &vb11, csq, snq, &r__);
#line 384 "zlags2.f"
	    }

#line 386 "zlags2.f"
	    *csu = snr;
#line 387 "zlags2.f"
	    d_cnjg(&z__2, &d1);
#line 387 "zlags2.f"
	    z__1.r = csr * z__2.r, z__1.i = csr * z__2.i;
#line 387 "zlags2.f"
	    snu->r = z__1.r, snu->i = z__1.i;
#line 388 "zlags2.f"
	    *csv = snl;
#line 389 "zlags2.f"
	    d_cnjg(&z__2, &d1);
#line 389 "zlags2.f"
	    z__1.r = csl * z__2.r, z__1.i = csl * z__2.i;
#line 389 "zlags2.f"
	    snv->r = z__1.r, snv->i = z__1.i;

#line 391 "zlags2.f"
	}

#line 393 "zlags2.f"
    }

#line 395 "zlags2.f"
    return 0;

/*     End of ZLAGS2 */

} /* zlags2_ */

