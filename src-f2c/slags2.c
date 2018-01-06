#line 1 "slags2.f"
/* slags2.f -- translated by f2c (version 20100827).
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

#line 1 "slags2.f"
/* > \brief \b SLAGS2 computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B su
ch that the rows of the transformed A and B are parallel. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slags2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slags2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slags2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, */
/*                          SNV, CSQ, SNQ ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            UPPER */
/*       REAL               A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, SNQ, */
/*      $                   SNU, SNV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such */
/* > that if ( UPPER ) then */
/* > */
/* >           U**T *A*Q = U**T *( A1 A2 )*Q = ( x  0  ) */
/* >                             ( 0  A3 )     ( x  x  ) */
/* > and */
/* >           V**T*B*Q = V**T *( B1 B2 )*Q = ( x  0  ) */
/* >                            ( 0  B3 )     ( x  x  ) */
/* > */
/* > or if ( .NOT.UPPER ) then */
/* > */
/* >           U**T *A*Q = U**T *( A1 0  )*Q = ( x  x  ) */
/* >                             ( A2 A3 )     ( 0  x  ) */
/* > and */
/* >           V**T*B*Q = V**T*( B1 0  )*Q = ( x  x  ) */
/* >                           ( B2 B3 )     ( 0  x  ) */
/* > */
/* > The rows of the transformed A and B are parallel, where */
/* > */
/* >   U = (  CSU  SNU ), V = (  CSV SNV ), Q = (  CSQ   SNQ ) */
/* >       ( -SNU  CSU )      ( -SNV CSV )      ( -SNQ   CSQ ) */
/* > */
/* > Z**T denotes the transpose of Z. */
/* > */
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
/* >          A1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] A2 */
/* > \verbatim */
/* >          A2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] A3 */
/* > \verbatim */
/* >          A3 is REAL */
/* >          On entry, A1, A2 and A3 are elements of the input 2-by-2 */
/* >          upper (lower) triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] B1 */
/* > \verbatim */
/* >          B1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B2 */
/* > \verbatim */
/* >          B2 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] B3 */
/* > \verbatim */
/* >          B3 is REAL */
/* >          On entry, B1, B2 and B3 are elements of the input 2-by-2 */
/* >          upper (lower) triangular matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] CSU */
/* > \verbatim */
/* >          CSU is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNU */
/* > \verbatim */
/* >          SNU is REAL */
/* >          The desired orthogonal matrix U. */
/* > \endverbatim */
/* > */
/* > \param[out] CSV */
/* > \verbatim */
/* >          CSV is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNV */
/* > \verbatim */
/* >          SNV is REAL */
/* >          The desired orthogonal matrix V. */
/* > \endverbatim */
/* > */
/* > \param[out] CSQ */
/* > \verbatim */
/* >          CSQ is REAL */
/* > \endverbatim */
/* > */
/* > \param[out] SNQ */
/* > \verbatim */
/* >          SNQ is REAL */
/* >          The desired orthogonal matrix Q. */
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
/* Subroutine */ int slags2_(logical *upper, doublereal *a1, doublereal *a2, 
	doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3, 
	doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv, 
	doublereal *csq, doublereal *snq)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal a, b, c__, d__, r__, s1, s2, ua11, ua12, ua21, ua22, 
	    vb11, vb12, vb21, vb22, csl, csr, snl, snr, aua11, aua12, aua21, 
	    aua22, avb11, avb12, avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    extern /* Subroutine */ int slasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), slartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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
/*     .. Executable Statements .. */

#line 186 "slags2.f"
    if (*upper) {

/*        Input matrices A and B are upper triangular matrices */

/*        Form matrix C = A*adj(B) = ( a b ) */
/*                                   ( 0 d ) */

#line 193 "slags2.f"
	a = *a1 * *b3;
#line 194 "slags2.f"
	d__ = *a3 * *b1;
#line 195 "slags2.f"
	b = *a2 * *b1 - *a1 * *b2;

/*        The SVD of real 2-by-2 triangular C */

/*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T ) */

#line 202 "slags2.f"
	slasv2_(&a, &b, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

#line 204 "slags2.f"
	if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {

/*           Compute the (1,1) and (1,2) elements of U**T *A and V**T *B, */
/*           and (1,2) element of |U|**T *|A| and |V|**T *|B|. */

#line 210 "slags2.f"
	    ua11r = csl * *a1;
#line 211 "slags2.f"
	    ua12 = csl * *a2 + snl * *a3;

#line 213 "slags2.f"
	    vb11r = csr * *b1;
#line 214 "slags2.f"
	    vb12 = csr * *b2 + snr * *b3;

#line 216 "slags2.f"
	    aua12 = abs(csl) * abs(*a2) + abs(snl) * abs(*a3);
#line 217 "slags2.f"
	    avb12 = abs(csr) * abs(*b2) + abs(snr) * abs(*b3);

/*           zero (1,2) elements of U**T *A and V**T *B */

#line 221 "slags2.f"
	    if (abs(ua11r) + abs(ua12) != 0.) {
#line 222 "slags2.f"
		if (aua12 / (abs(ua11r) + abs(ua12)) <= avb12 / (abs(vb11r) + 
			abs(vb12))) {
#line 224 "slags2.f"
		    d__1 = -ua11r;
#line 224 "slags2.f"
		    slartg_(&d__1, &ua12, csq, snq, &r__);
#line 225 "slags2.f"
		} else {
#line 226 "slags2.f"
		    d__1 = -vb11r;
#line 226 "slags2.f"
		    slartg_(&d__1, &vb12, csq, snq, &r__);
#line 227 "slags2.f"
		}
#line 228 "slags2.f"
	    } else {
#line 229 "slags2.f"
		d__1 = -vb11r;
#line 229 "slags2.f"
		slartg_(&d__1, &vb12, csq, snq, &r__);
#line 230 "slags2.f"
	    }

#line 232 "slags2.f"
	    *csu = csl;
#line 233 "slags2.f"
	    *snu = -snl;
#line 234 "slags2.f"
	    *csv = csr;
#line 235 "slags2.f"
	    *snv = -snr;

#line 237 "slags2.f"
	} else {

/*           Compute the (2,1) and (2,2) elements of U**T *A and V**T *B, */
/*           and (2,2) element of |U|**T *|A| and |V|**T *|B|. */

#line 242 "slags2.f"
	    ua21 = -snl * *a1;
#line 243 "slags2.f"
	    ua22 = -snl * *a2 + csl * *a3;

#line 245 "slags2.f"
	    vb21 = -snr * *b1;
#line 246 "slags2.f"
	    vb22 = -snr * *b2 + csr * *b3;

#line 248 "slags2.f"
	    aua22 = abs(snl) * abs(*a2) + abs(csl) * abs(*a3);
#line 249 "slags2.f"
	    avb22 = abs(snr) * abs(*b2) + abs(csr) * abs(*b3);

/*           zero (2,2) elements of U**T*A and V**T*B, and then swap. */

#line 253 "slags2.f"
	    if (abs(ua21) + abs(ua22) != 0.) {
#line 254 "slags2.f"
		if (aua22 / (abs(ua21) + abs(ua22)) <= avb22 / (abs(vb21) + 
			abs(vb22))) {
#line 256 "slags2.f"
		    d__1 = -ua21;
#line 256 "slags2.f"
		    slartg_(&d__1, &ua22, csq, snq, &r__);
#line 257 "slags2.f"
		} else {
#line 258 "slags2.f"
		    d__1 = -vb21;
#line 258 "slags2.f"
		    slartg_(&d__1, &vb22, csq, snq, &r__);
#line 259 "slags2.f"
		}
#line 260 "slags2.f"
	    } else {
#line 261 "slags2.f"
		d__1 = -vb21;
#line 261 "slags2.f"
		slartg_(&d__1, &vb22, csq, snq, &r__);
#line 262 "slags2.f"
	    }

#line 264 "slags2.f"
	    *csu = snl;
#line 265 "slags2.f"
	    *snu = csl;
#line 266 "slags2.f"
	    *csv = snr;
#line 267 "slags2.f"
	    *snv = csr;

#line 269 "slags2.f"
	}

#line 271 "slags2.f"
    } else {

/*        Input matrices A and B are lower triangular matrices */

/*        Form matrix C = A*adj(B) = ( a 0 ) */
/*                                   ( c d ) */

#line 278 "slags2.f"
	a = *a1 * *b3;
#line 279 "slags2.f"
	d__ = *a3 * *b1;
#line 280 "slags2.f"
	c__ = *a2 * *b3 - *a3 * *b2;

/*        The SVD of real 2-by-2 triangular C */

/*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 ) */
/*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T ) */

#line 287 "slags2.f"
	slasv2_(&a, &c__, &d__, &s1, &s2, &snr, &csr, &snl, &csl);

#line 289 "slags2.f"
	if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {

/*           Compute the (2,1) and (2,2) elements of U**T *A and V**T *B, */
/*           and (2,1) element of |U|**T *|A| and |V|**T *|B|. */

#line 295 "slags2.f"
	    ua21 = -snr * *a1 + csr * *a2;
#line 296 "slags2.f"
	    ua22r = csr * *a3;

#line 298 "slags2.f"
	    vb21 = -snl * *b1 + csl * *b2;
#line 299 "slags2.f"
	    vb22r = csl * *b3;

#line 301 "slags2.f"
	    aua21 = abs(snr) * abs(*a1) + abs(csr) * abs(*a2);
#line 302 "slags2.f"
	    avb21 = abs(snl) * abs(*b1) + abs(csl) * abs(*b2);

/*           zero (2,1) elements of U**T *A and V**T *B. */

#line 306 "slags2.f"
	    if (abs(ua21) + abs(ua22r) != 0.) {
#line 307 "slags2.f"
		if (aua21 / (abs(ua21) + abs(ua22r)) <= avb21 / (abs(vb21) + 
			abs(vb22r))) {
#line 309 "slags2.f"
		    slartg_(&ua22r, &ua21, csq, snq, &r__);
#line 310 "slags2.f"
		} else {
#line 311 "slags2.f"
		    slartg_(&vb22r, &vb21, csq, snq, &r__);
#line 312 "slags2.f"
		}
#line 313 "slags2.f"
	    } else {
#line 314 "slags2.f"
		slartg_(&vb22r, &vb21, csq, snq, &r__);
#line 315 "slags2.f"
	    }

#line 317 "slags2.f"
	    *csu = csr;
#line 318 "slags2.f"
	    *snu = -snr;
#line 319 "slags2.f"
	    *csv = csl;
#line 320 "slags2.f"
	    *snv = -snl;

#line 322 "slags2.f"
	} else {

/*           Compute the (1,1) and (1,2) elements of U**T *A and V**T *B, */
/*           and (1,1) element of |U|**T *|A| and |V|**T *|B|. */

#line 327 "slags2.f"
	    ua11 = csr * *a1 + snr * *a2;
#line 328 "slags2.f"
	    ua12 = snr * *a3;

#line 330 "slags2.f"
	    vb11 = csl * *b1 + snl * *b2;
#line 331 "slags2.f"
	    vb12 = snl * *b3;

#line 333 "slags2.f"
	    aua11 = abs(csr) * abs(*a1) + abs(snr) * abs(*a2);
#line 334 "slags2.f"
	    avb11 = abs(csl) * abs(*b1) + abs(snl) * abs(*b2);

/*           zero (1,1) elements of U**T*A and V**T*B, and then swap. */

#line 338 "slags2.f"
	    if (abs(ua11) + abs(ua12) != 0.) {
#line 339 "slags2.f"
		if (aua11 / (abs(ua11) + abs(ua12)) <= avb11 / (abs(vb11) + 
			abs(vb12))) {
#line 341 "slags2.f"
		    slartg_(&ua12, &ua11, csq, snq, &r__);
#line 342 "slags2.f"
		} else {
#line 343 "slags2.f"
		    slartg_(&vb12, &vb11, csq, snq, &r__);
#line 344 "slags2.f"
		}
#line 345 "slags2.f"
	    } else {
#line 346 "slags2.f"
		slartg_(&vb12, &vb11, csq, snq, &r__);
#line 347 "slags2.f"
	    }

#line 349 "slags2.f"
	    *csu = snr;
#line 350 "slags2.f"
	    *snu = csr;
#line 351 "slags2.f"
	    *csv = snl;
#line 352 "slags2.f"
	    *snv = csl;

#line 354 "slags2.f"
	}

#line 356 "slags2.f"
    }

#line 358 "slags2.f"
    return 0;

/*     End of SLAGS2 */

} /* slags2_ */

