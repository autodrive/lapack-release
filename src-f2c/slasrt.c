#line 1 "slasrt.f"
/* slasrt.f -- translated by f2c (version 20100827).
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

#line 1 "slasrt.f"
/* > \brief \b SLASRT sorts numbers in increasing or decreasing order. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasrt.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasrt.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasrt.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASRT( ID, N, D, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          ID */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Sort the numbers in D in increasing order (if ID = 'I') or */
/* > in decreasing order (if ID = 'D' ). */
/* > */
/* > Use Quick Sort, reverting to Insertion sort on arrays of */
/* > size <= 20. Dimension of STACK limits N to about 2**32. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ID */
/* > \verbatim */
/* >          ID is CHARACTER*1 */
/* >          = 'I': sort D in increasing order; */
/* >          = 'D': sort D in decreasing order. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The length of the array D. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          On entry, the array to be sorted. */
/* >          On exit, D has been sorted into increasing order */
/* >          (D(1) <= ... <= D(N) ) or into decreasing order */
/* >          (D(1) >= ... >= D(N) ), depending on ID. */
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

/* > \date June 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int slasrt_(char *id, integer *n, doublereal *d__, integer *
	info, ftnlen id_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal d1, d2, d3;
    static integer dir;
    static doublereal tmp;
    static integer endd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer stack[64]	/* was [2][32] */;
    static doublereal dmnmx;
    static integer start;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer stkpnt;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 128 "slasrt.f"
    /* Parameter adjustments */
#line 128 "slasrt.f"
    --d__;
#line 128 "slasrt.f"

#line 128 "slasrt.f"
    /* Function Body */
#line 128 "slasrt.f"
    *info = 0;
#line 129 "slasrt.f"
    dir = -1;
#line 130 "slasrt.f"
    if (lsame_(id, "D", (ftnlen)1, (ftnlen)1)) {
#line 131 "slasrt.f"
	dir = 0;
#line 132 "slasrt.f"
    } else if (lsame_(id, "I", (ftnlen)1, (ftnlen)1)) {
#line 133 "slasrt.f"
	dir = 1;
#line 134 "slasrt.f"
    }
#line 135 "slasrt.f"
    if (dir == -1) {
#line 136 "slasrt.f"
	*info = -1;
#line 137 "slasrt.f"
    } else if (*n < 0) {
#line 138 "slasrt.f"
	*info = -2;
#line 139 "slasrt.f"
    }
#line 140 "slasrt.f"
    if (*info != 0) {
#line 141 "slasrt.f"
	i__1 = -(*info);
#line 141 "slasrt.f"
	xerbla_("SLASRT", &i__1, (ftnlen)6);
#line 142 "slasrt.f"
	return 0;
#line 143 "slasrt.f"
    }

/*     Quick return if possible */

#line 147 "slasrt.f"
    if (*n <= 1) {
#line 147 "slasrt.f"
	return 0;
#line 147 "slasrt.f"
    }

#line 150 "slasrt.f"
    stkpnt = 1;
#line 151 "slasrt.f"
    stack[0] = 1;
#line 152 "slasrt.f"
    stack[1] = *n;
#line 153 "slasrt.f"
L10:
#line 154 "slasrt.f"
    start = stack[(stkpnt << 1) - 2];
#line 155 "slasrt.f"
    endd = stack[(stkpnt << 1) - 1];
#line 156 "slasrt.f"
    --stkpnt;
#line 157 "slasrt.f"
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

#line 161 "slasrt.f"
	if (dir == 0) {

/*           Sort into decreasing order */

#line 165 "slasrt.f"
	    i__1 = endd;
#line 165 "slasrt.f"
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
#line 166 "slasrt.f"
		i__2 = start + 1;
#line 166 "slasrt.f"
		for (j = i__; j >= i__2; --j) {
#line 167 "slasrt.f"
		    if (d__[j] > d__[j - 1]) {
#line 168 "slasrt.f"
			dmnmx = d__[j];
#line 169 "slasrt.f"
			d__[j] = d__[j - 1];
#line 170 "slasrt.f"
			d__[j - 1] = dmnmx;
#line 171 "slasrt.f"
		    } else {
#line 172 "slasrt.f"
			goto L30;
#line 173 "slasrt.f"
		    }
#line 174 "slasrt.f"
/* L20: */
#line 174 "slasrt.f"
		}
#line 175 "slasrt.f"
L30:
#line 175 "slasrt.f"
		;
#line 175 "slasrt.f"
	    }

#line 177 "slasrt.f"
	} else {

/*           Sort into increasing order */

#line 181 "slasrt.f"
	    i__1 = endd;
#line 181 "slasrt.f"
	    for (i__ = start + 1; i__ <= i__1; ++i__) {
#line 182 "slasrt.f"
		i__2 = start + 1;
#line 182 "slasrt.f"
		for (j = i__; j >= i__2; --j) {
#line 183 "slasrt.f"
		    if (d__[j] < d__[j - 1]) {
#line 184 "slasrt.f"
			dmnmx = d__[j];
#line 185 "slasrt.f"
			d__[j] = d__[j - 1];
#line 186 "slasrt.f"
			d__[j - 1] = dmnmx;
#line 187 "slasrt.f"
		    } else {
#line 188 "slasrt.f"
			goto L50;
#line 189 "slasrt.f"
		    }
#line 190 "slasrt.f"
/* L40: */
#line 190 "slasrt.f"
		}
#line 191 "slasrt.f"
L50:
#line 191 "slasrt.f"
		;
#line 191 "slasrt.f"
	    }

#line 193 "slasrt.f"
	}

#line 195 "slasrt.f"
    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first */

/*        Choose partition entry as median of 3 */

#line 201 "slasrt.f"
	d1 = d__[start];
#line 202 "slasrt.f"
	d2 = d__[endd];
#line 203 "slasrt.f"
	i__ = (start + endd) / 2;
#line 204 "slasrt.f"
	d3 = d__[i__];
#line 205 "slasrt.f"
	if (d1 < d2) {
#line 206 "slasrt.f"
	    if (d3 < d1) {
#line 207 "slasrt.f"
		dmnmx = d1;
#line 208 "slasrt.f"
	    } else if (d3 < d2) {
#line 209 "slasrt.f"
		dmnmx = d3;
#line 210 "slasrt.f"
	    } else {
#line 211 "slasrt.f"
		dmnmx = d2;
#line 212 "slasrt.f"
	    }
#line 213 "slasrt.f"
	} else {
#line 214 "slasrt.f"
	    if (d3 < d2) {
#line 215 "slasrt.f"
		dmnmx = d2;
#line 216 "slasrt.f"
	    } else if (d3 < d1) {
#line 217 "slasrt.f"
		dmnmx = d3;
#line 218 "slasrt.f"
	    } else {
#line 219 "slasrt.f"
		dmnmx = d1;
#line 220 "slasrt.f"
	    }
#line 221 "slasrt.f"
	}

#line 223 "slasrt.f"
	if (dir == 0) {

/*           Sort into decreasing order */

#line 227 "slasrt.f"
	    i__ = start - 1;
#line 228 "slasrt.f"
	    j = endd + 1;
#line 229 "slasrt.f"
L60:
#line 230 "slasrt.f"
L70:
#line 231 "slasrt.f"
	    --j;
#line 232 "slasrt.f"
	    if (d__[j] < dmnmx) {
#line 232 "slasrt.f"
		goto L70;
#line 232 "slasrt.f"
	    }
#line 234 "slasrt.f"
L80:
#line 235 "slasrt.f"
	    ++i__;
#line 236 "slasrt.f"
	    if (d__[i__] > dmnmx) {
#line 236 "slasrt.f"
		goto L80;
#line 236 "slasrt.f"
	    }
#line 238 "slasrt.f"
	    if (i__ < j) {
#line 239 "slasrt.f"
		tmp = d__[i__];
#line 240 "slasrt.f"
		d__[i__] = d__[j];
#line 241 "slasrt.f"
		d__[j] = tmp;
#line 242 "slasrt.f"
		goto L60;
#line 243 "slasrt.f"
	    }
#line 244 "slasrt.f"
	    if (j - start > endd - j - 1) {
#line 245 "slasrt.f"
		++stkpnt;
#line 246 "slasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 247 "slasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 248 "slasrt.f"
		++stkpnt;
#line 249 "slasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 250 "slasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 251 "slasrt.f"
	    } else {
#line 252 "slasrt.f"
		++stkpnt;
#line 253 "slasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 254 "slasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 255 "slasrt.f"
		++stkpnt;
#line 256 "slasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 257 "slasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 258 "slasrt.f"
	    }
#line 259 "slasrt.f"
	} else {

/*           Sort into increasing order */

#line 263 "slasrt.f"
	    i__ = start - 1;
#line 264 "slasrt.f"
	    j = endd + 1;
#line 265 "slasrt.f"
L90:
#line 266 "slasrt.f"
L100:
#line 267 "slasrt.f"
	    --j;
#line 268 "slasrt.f"
	    if (d__[j] > dmnmx) {
#line 268 "slasrt.f"
		goto L100;
#line 268 "slasrt.f"
	    }
#line 270 "slasrt.f"
L110:
#line 271 "slasrt.f"
	    ++i__;
#line 272 "slasrt.f"
	    if (d__[i__] < dmnmx) {
#line 272 "slasrt.f"
		goto L110;
#line 272 "slasrt.f"
	    }
#line 274 "slasrt.f"
	    if (i__ < j) {
#line 275 "slasrt.f"
		tmp = d__[i__];
#line 276 "slasrt.f"
		d__[i__] = d__[j];
#line 277 "slasrt.f"
		d__[j] = tmp;
#line 278 "slasrt.f"
		goto L90;
#line 279 "slasrt.f"
	    }
#line 280 "slasrt.f"
	    if (j - start > endd - j - 1) {
#line 281 "slasrt.f"
		++stkpnt;
#line 282 "slasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 283 "slasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 284 "slasrt.f"
		++stkpnt;
#line 285 "slasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 286 "slasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 287 "slasrt.f"
	    } else {
#line 288 "slasrt.f"
		++stkpnt;
#line 289 "slasrt.f"
		stack[(stkpnt << 1) - 2] = j + 1;
#line 290 "slasrt.f"
		stack[(stkpnt << 1) - 1] = endd;
#line 291 "slasrt.f"
		++stkpnt;
#line 292 "slasrt.f"
		stack[(stkpnt << 1) - 2] = start;
#line 293 "slasrt.f"
		stack[(stkpnt << 1) - 1] = j;
#line 294 "slasrt.f"
	    }
#line 295 "slasrt.f"
	}
#line 296 "slasrt.f"
    }
#line 297 "slasrt.f"
    if (stkpnt > 0) {
#line 297 "slasrt.f"
	goto L10;
#line 297 "slasrt.f"
    }
#line 299 "slasrt.f"
    return 0;

/*     End of SLASRT */

} /* slasrt_ */

