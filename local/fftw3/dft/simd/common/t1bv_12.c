/*
 * Copyright (c) 2003, 2007-11 Matteo Frigo
 * Copyright (c) 2003, 2007-11 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Wed Jul 27 06:15:34 EDT 2011 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c.native -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 12 -name t1bv_12 -include t1b.h -sign 1 */

/*
 * This function contains 59 FP additions, 42 FP multiplications,
 * (or, 41 additions, 24 multiplications, 18 fused multiply/add),
 * 41 stack variables, 2 constants, and 24 memory accesses
 */
#include "t1b.h"

static void t1bv_12(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  R *x;
	  x = ii;
	  for (m = mb, W = W + (mb * ((TWVL / VL) * 22)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 22), MAKE_VOLATILE_STRIDE(rs)) {
	       V TI, Ti, TA, T7, Tm, TE, Tw, Tk, Tf, TB, TU, TM;
	       {
		    V T9, TK, Tj, TL, Te;
		    {
			 V T1, T4, T2, Tp, Tt, Tr;
			 T1 = LD(&(x[0]), ms, &(x[0]));
			 T4 = LD(&(x[WS(rs, 8)]), ms, &(x[0]));
			 T2 = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
			 Tp = LD(&(x[WS(rs, 9)]), ms, &(x[WS(rs, 1)]));
			 Tt = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
			 Tr = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
			 {
			      V T5, T3, Tq, Tu, Ts, Td, Tb, T8, Tc, Ta;
			      T8 = LD(&(x[WS(rs, 6)]), ms, &(x[0]));
			      Tc = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
			      Ta = LD(&(x[WS(rs, 10)]), ms, &(x[0]));
			      T5 = BYTW(&(W[TWVL * 14]), T4);
			      T3 = BYTW(&(W[TWVL * 6]), T2);
			      Tq = BYTW(&(W[TWVL * 16]), Tp);
			      Tu = BYTW(&(W[TWVL * 8]), Tt);
			      Ts = BYTW(&(W[0]), Tr);
			      T9 = BYTW(&(W[TWVL * 10]), T8);
			      Td = BYTW(&(W[TWVL * 2]), Tc);
			      Tb = BYTW(&(W[TWVL * 18]), Ta);
			      {
				   V Th, T6, Tl, Tv;
				   Th = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
				   TK = VSUB(T3, T5);
				   T6 = VADD(T3, T5);
				   Tl = LD(&(x[WS(rs, 11)]), ms, &(x[WS(rs, 1)]));
				   Tv = VADD(Ts, Tu);
				   TI = VSUB(Tu, Ts);
				   Tj = LD(&(x[WS(rs, 7)]), ms, &(x[WS(rs, 1)]));
				   TL = VSUB(Tb, Td);
				   Te = VADD(Tb, Td);
				   Ti = BYTW(&(W[TWVL * 4]), Th);
				   TA = VFNMS(LDK(KP500000000), T6, T1);
				   T7 = VADD(T1, T6);
				   Tm = BYTW(&(W[TWVL * 20]), Tl);
				   TE = VFNMS(LDK(KP500000000), Tv, Tq);
				   Tw = VADD(Tq, Tv);
			      }
			 }
		    }
		    Tk = BYTW(&(W[TWVL * 12]), Tj);
		    Tf = VADD(T9, Te);
		    TB = VFNMS(LDK(KP500000000), Te, T9);
		    TU = VSUB(TK, TL);
		    TM = VADD(TK, TL);
	       }
	       {
		    V Tn, TH, TC, TQ, Ty, Tg;
		    Tn = VADD(Tk, Tm);
		    TH = VSUB(Tk, Tm);
		    TC = VADD(TA, TB);
		    TQ = VSUB(TA, TB);
		    Ty = VADD(T7, Tf);
		    Tg = VSUB(T7, Tf);
		    {
			 V To, TD, TJ, TR;
			 To = VADD(Ti, Tn);
			 TD = VFNMS(LDK(KP500000000), Tn, Ti);
			 TJ = VSUB(TH, TI);
			 TR = VADD(TH, TI);
			 {
			      V TP, TN, TW, TS, TO, TG, TX, TV;
			      {
				   V Tz, Tx, TF, TT;
				   Tz = VADD(To, Tw);
				   Tx = VSUB(To, Tw);
				   TF = VADD(TD, TE);
				   TT = VSUB(TD, TE);
				   TP = VMUL(LDK(KP866025403), VADD(TM, TJ));
				   TN = VMUL(LDK(KP866025403), VSUB(TJ, TM));
				   TW = VFMA(LDK(KP866025403), TR, TQ);
				   TS = VFNMS(LDK(KP866025403), TR, TQ);
				   ST(&(x[WS(rs, 6)]), VSUB(Ty, Tz), ms, &(x[0]));
				   ST(&(x[0]), VADD(Ty, Tz), ms, &(x[0]));
				   ST(&(x[WS(rs, 9)]), VFMAI(Tx, Tg), ms, &(x[WS(rs, 1)]));
				   ST(&(x[WS(rs, 3)]), VFNMSI(Tx, Tg), ms, &(x[WS(rs, 1)]));
				   TO = VADD(TC, TF);
				   TG = VSUB(TC, TF);
				   TX = VFNMS(LDK(KP866025403), TU, TT);
				   TV = VFMA(LDK(KP866025403), TU, TT);
			      }
			      ST(&(x[WS(rs, 8)]), VFNMSI(TP, TO), ms, &(x[0]));
			      ST(&(x[WS(rs, 4)]), VFMAI(TP, TO), ms, &(x[0]));
			      ST(&(x[WS(rs, 2)]), VFMAI(TN, TG), ms, &(x[0]));
			      ST(&(x[WS(rs, 10)]), VFNMSI(TN, TG), ms, &(x[0]));
			      ST(&(x[WS(rs, 5)]), VFMAI(TX, TW), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 7)]), VFNMSI(TX, TW), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 11)]), VFNMSI(TV, TS), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 1)]), VFMAI(TV, TS), ms, &(x[WS(rs, 1)]));
			 }
		    }
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 2),
     VTW(0, 3),
     VTW(0, 4),
     VTW(0, 5),
     VTW(0, 6),
     VTW(0, 7),
     VTW(0, 8),
     VTW(0, 9),
     VTW(0, 10),
     VTW(0, 11),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 12, XSIMD_STRING("t1bv_12"), twinstr, &GENUS, {41, 24, 18, 0}, 0, 0, 0 };

void XSIMD(codelet_t1bv_12) (planner *p) {
     X(kdft_dit_register) (p, t1bv_12, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c.native -simd -compact -variables 4 -pipeline-latency 8 -n 12 -name t1bv_12 -include t1b.h -sign 1 */

/*
 * This function contains 59 FP additions, 30 FP multiplications,
 * (or, 55 additions, 26 multiplications, 4 fused multiply/add),
 * 28 stack variables, 2 constants, and 24 memory accesses
 */
#include "t1b.h"

static void t1bv_12(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DVK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  R *x;
	  x = ii;
	  for (m = mb, W = W + (mb * ((TWVL / VL) * 22)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 22), MAKE_VOLATILE_STRIDE(rs)) {
	       V T1, Tt, T6, T7, TB, Tq, TC, TD, T9, Tu, Te, Tf, Tx, Tl, Ty;
	       V Tz;
	       {
		    V T5, T3, T4, T2;
		    T1 = LD(&(x[0]), ms, &(x[0]));
		    T4 = LD(&(x[WS(rs, 8)]), ms, &(x[0]));
		    T5 = BYTW(&(W[TWVL * 14]), T4);
		    T2 = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
		    T3 = BYTW(&(W[TWVL * 6]), T2);
		    Tt = VSUB(T3, T5);
		    T6 = VADD(T3, T5);
		    T7 = VFNMS(LDK(KP500000000), T6, T1);
	       }
	       {
		    V Tn, Tp, Tm, TA, To;
		    Tm = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
		    Tn = BYTW(&(W[0]), Tm);
		    TA = LD(&(x[WS(rs, 9)]), ms, &(x[WS(rs, 1)]));
		    TB = BYTW(&(W[TWVL * 16]), TA);
		    To = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
		    Tp = BYTW(&(W[TWVL * 8]), To);
		    Tq = VSUB(Tn, Tp);
		    TC = VADD(Tn, Tp);
		    TD = VFNMS(LDK(KP500000000), TC, TB);
	       }
	       {
		    V Td, Tb, T8, Tc, Ta;
		    T8 = LD(&(x[WS(rs, 6)]), ms, &(x[0]));
		    T9 = BYTW(&(W[TWVL * 10]), T8);
		    Tc = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
		    Td = BYTW(&(W[TWVL * 2]), Tc);
		    Ta = LD(&(x[WS(rs, 10)]), ms, &(x[0]));
		    Tb = BYTW(&(W[TWVL * 18]), Ta);
		    Tu = VSUB(Tb, Td);
		    Te = VADD(Tb, Td);
		    Tf = VFNMS(LDK(KP500000000), Te, T9);
	       }
	       {
		    V Ti, Tk, Th, Tw, Tj;
		    Th = LD(&(x[WS(rs, 7)]), ms, &(x[WS(rs, 1)]));
		    Ti = BYTW(&(W[TWVL * 12]), Th);
		    Tw = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
		    Tx = BYTW(&(W[TWVL * 4]), Tw);
		    Tj = LD(&(x[WS(rs, 11)]), ms, &(x[WS(rs, 1)]));
		    Tk = BYTW(&(W[TWVL * 20]), Tj);
		    Tl = VSUB(Ti, Tk);
		    Ty = VADD(Ti, Tk);
		    Tz = VFNMS(LDK(KP500000000), Ty, Tx);
	       }
	       {
		    V Ts, TG, TF, TH;
		    {
			 V Tg, Tr, Tv, TE;
			 Tg = VSUB(T7, Tf);
			 Tr = VMUL(LDK(KP866025403), VSUB(Tl, Tq));
			 Ts = VSUB(Tg, Tr);
			 TG = VADD(Tg, Tr);
			 Tv = VMUL(LDK(KP866025403), VSUB(Tt, Tu));
			 TE = VSUB(Tz, TD);
			 TF = VBYI(VADD(Tv, TE));
			 TH = VBYI(VSUB(TE, Tv));
		    }
		    ST(&(x[WS(rs, 11)]), VSUB(Ts, TF), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 5)]), VADD(TG, TH), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 1)]), VADD(Ts, TF), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 7)]), VSUB(TG, TH), ms, &(x[WS(rs, 1)]));
	       }
	       {
		    V TS, TW, TV, TX;
		    {
			 V TQ, TR, TT, TU;
			 TQ = VADD(T1, T6);
			 TR = VADD(T9, Te);
			 TS = VSUB(TQ, TR);
			 TW = VADD(TQ, TR);
			 TT = VADD(Tx, Ty);
			 TU = VADD(TB, TC);
			 TV = VBYI(VSUB(TT, TU));
			 TX = VADD(TT, TU);
		    }
		    ST(&(x[WS(rs, 3)]), VSUB(TS, TV), ms, &(x[WS(rs, 1)]));
		    ST(&(x[0]), VADD(TW, TX), ms, &(x[0]));
		    ST(&(x[WS(rs, 9)]), VADD(TS, TV), ms, &(x[WS(rs, 1)]));
		    ST(&(x[WS(rs, 6)]), VSUB(TW, TX), ms, &(x[0]));
	       }
	       {
		    V TK, TO, TN, TP;
		    {
			 V TI, TJ, TL, TM;
			 TI = VADD(Tl, Tq);
			 TJ = VADD(Tt, Tu);
			 TK = VBYI(VMUL(LDK(KP866025403), VSUB(TI, TJ)));
			 TO = VBYI(VMUL(LDK(KP866025403), VADD(TJ, TI)));
			 TL = VADD(T7, Tf);
			 TM = VADD(Tz, TD);
			 TN = VSUB(TL, TM);
			 TP = VADD(TL, TM);
		    }
		    ST(&(x[WS(rs, 2)]), VADD(TK, TN), ms, &(x[0]));
		    ST(&(x[WS(rs, 8)]), VSUB(TP, TO), ms, &(x[0]));
		    ST(&(x[WS(rs, 10)]), VSUB(TN, TK), ms, &(x[0]));
		    ST(&(x[WS(rs, 4)]), VADD(TO, TP), ms, &(x[0]));
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 2),
     VTW(0, 3),
     VTW(0, 4),
     VTW(0, 5),
     VTW(0, 6),
     VTW(0, 7),
     VTW(0, 8),
     VTW(0, 9),
     VTW(0, 10),
     VTW(0, 11),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 12, XSIMD_STRING("t1bv_12"), twinstr, &GENUS, {55, 26, 4, 0}, 0, 0, 0 };

void XSIMD(codelet_t1bv_12) (planner *p) {
     X(kdft_dit_register) (p, t1bv_12, &desc);
}
#endif				/* HAVE_FMA */
