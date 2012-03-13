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
/* Generated on Wed Jul 27 06:15:25 EDT 2011 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c.native -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -twiddle-log3 -precompute-twiddles -no-generate-bytw -n 32 -name t3fv_32 -include t3f.h */

/*
 * This function contains 244 FP additions, 214 FP multiplications,
 * (or, 146 additions, 116 multiplications, 98 fused multiply/add),
 * 118 stack variables, 7 constants, and 64 memory accesses
 */
#include "t3f.h"

static void t3fv_32(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DVK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     {
	  INT m;
	  R *x;
	  x = ri;
	  for (m = mb, W = W + (mb * ((TWVL / VL) * 8)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 8), MAKE_VOLATILE_STRIDE(rs)) {
	       V T2B, T2A, T2u, T2x, T2r, T2F, T2L, T2P;
	       {
		    V T2, T5, T3, T7;
		    T2 = LDW(&(W[0]));
		    T5 = LDW(&(W[TWVL * 4]));
		    T3 = LDW(&(W[TWVL * 2]));
		    T7 = LDW(&(W[TWVL * 6]));
		    {
			 V T24, Tb, T3x, T2T, T3K, T2W, T25, Tr, T3z, T3g, T28, TX, T3y, T3j, T27;
			 V TG, T37, T3F, T3G, T3a, T2Y, T15, T1p, T2Z, T2w, T1V, T2v, T1N, T32, T1h;
			 V T17, T1a;
			 {
			      V T1, Tz, TT, T4, TC, Tv, T12, T1D, T1w, T18, T1t, T1O, TK, TP, T1c;
			      V T1m, Tf, T6, Te, TL, TQ, T2S, Tp, TU, Ti, Ta, TM, TR, Tm, TJ;
			      V T22, T9, T1Z;
			      T1 = LD(&(x[0]), ms, &(x[0]));
			      T22 = LD(&(x[WS(rs, 24)]), ms, &(x[0]));
			      T9 = LD(&(x[WS(rs, 16)]), ms, &(x[0]));
			      T1Z = LD(&(x[WS(rs, 8)]), ms, &(x[0]));
			      {
				   V Tn, TH, Tk, To, Th, Tg, T8, Tl, T20, T23, TI;
				   {
					V Td, T1C, Tc, T21;
					Td = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
					Tz = VZMUL(T2, T5);
					T1C = VZMULJ(T2, T5);
					Tn = VZMUL(T3, T5);
					TT = VZMULJ(T3, T5);
					Tc = VZMUL(T2, T3);
					T4 = VZMULJ(T2, T3);
					TH = VZMUL(T3, T7);
					T21 = VZMULJ(T3, T7);
					Tk = VZMUL(T2, T7);
					TC = VZMULJ(T2, T7);
					Tv = VZMULJ(T5, T7);
					T12 = VZMULJ(Tz, T7);
					T20 = VZMULJ(T1C, T1Z);
					T1D = VZMULJ(T1C, T7);
					T1w = VZMULJ(Tn, T7);
					T18 = VZMULJ(TT, T7);
					T1t = VZMUL(Tc, T7);
					T1O = VZMULJ(Tc, T7);
					TK = VZMUL(Tc, T5);
					TP = VZMULJ(Tc, T5);
					T1c = VZMUL(T4, T7);
					T1m = VZMULJ(T4, T7);
					Tf = VZMULJ(T4, T5);
					T6 = VZMUL(T4, T5);
					T23 = VZMULJ(T21, T22);
					Te = VZMULJ(Tc, Td);
				   }
				   TL = VZMULJ(TK, T7);
				   TQ = VZMULJ(TP, T7);
				   To = LD(&(x[WS(rs, 12)]), ms, &(x[0]));
				   Th = LD(&(x[WS(rs, 20)]), ms, &(x[0]));
				   Tg = VZMULJ(Tf, T7);
				   T8 = VZMULJ(T6, T7);
				   T2S = VADD(T20, T23);
				   T24 = VSUB(T20, T23);
				   Tl = LD(&(x[WS(rs, 28)]), ms, &(x[0]));
				   TI = LD(&(x[WS(rs, 30)]), ms, &(x[0]));
				   Tp = VZMULJ(Tn, To);
				   TU = LD(&(x[WS(rs, 6)]), ms, &(x[0]));
				   Ti = VZMULJ(Tg, Th);
				   Ta = VZMULJ(T8, T9);
				   TM = LD(&(x[WS(rs, 14)]), ms, &(x[0]));
				   TR = LD(&(x[WS(rs, 22)]), ms, &(x[0]));
				   Tm = VZMULJ(Tk, Tl);
				   TJ = VZMULJ(TH, TI);
			      }
			      {
				   V Tu, TE, Tw, TA;
				   {
					V T3e, TO, T3f, TW;
					{
					     V TV, T2U, Tj, T2R, TN, TS, T2V, Tq, Tt, TD;
					     Tt = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
					     TV = VZMULJ(TT, TU);
					     T2U = VADD(Te, Ti);
					     Tj = VSUB(Te, Ti);
					     T2R = VADD(T1, Ta);
					     Tb = VSUB(T1, Ta);
					     TN = VZMULJ(TL, TM);
					     TS = VZMULJ(TQ, TR);
					     T2V = VADD(Tm, Tp);
					     Tq = VSUB(Tm, Tp);
					     Tu = VZMULJ(T4, Tt);
					     TD = LD(&(x[WS(rs, 26)]), ms, &(x[0]));
					     T3x = VSUB(T2R, T2S);
					     T2T = VADD(T2R, T2S);
					     T3e = VADD(TJ, TN);
					     TO = VSUB(TJ, TN);
					     T3f = VADD(TV, TS);
					     TW = VSUB(TS, TV);
					     T3K = VSUB(T2V, T2U);
					     T2W = VADD(T2U, T2V);
					     T25 = VSUB(Tq, Tj);
					     Tr = VADD(Tj, Tq);
					     TE = VZMULJ(TC, TD);
					}
					Tw = LD(&(x[WS(rs, 18)]), ms, &(x[0]));
					T3z = VSUB(T3e, T3f);
					T3g = VADD(T3e, T3f);
					T28 = VFMA(LDK(KP414213562), TO, TW);
					TX = VFNMS(LDK(KP414213562), TW, TO);
					TA = LD(&(x[WS(rs, 10)]), ms, &(x[0]));
				   }
				   {
					V T35, T1z, T1T, T36, T39, T1L, T1B, T1F;
					{
					     V T1v, T1y, Ty, T3h, T1S, T1Q, T1I, T3i, TF, T1K, T1A, T1E;
					     {
						  V T1u, T1x, Tx, T1R;
						  T1u = LD(&(x[WS(rs, 31)]), ms, &(x[WS(rs, 1)]));
						  T1x = LD(&(x[WS(rs, 15)]), ms, &(x[WS(rs, 1)]));
						  Tx = VZMULJ(Tv, Tw);
						  T1R = LD(&(x[WS(rs, 7)]), ms, &(x[WS(rs, 1)]));
						  {
						       V T1P, T1H, T1J, TB;
						       T1P = LD(&(x[WS(rs, 23)]), ms, &(x[WS(rs, 1)]));
						       T1H = LD(&(x[WS(rs, 27)]), ms, &(x[WS(rs, 1)]));
						       T1J = LD(&(x[WS(rs, 11)]), ms, &(x[WS(rs, 1)]));
						       TB = VZMULJ(Tz, TA);
						       T1v = VZMULJ(T1t, T1u);
						       T1y = VZMULJ(T1w, T1x);
						       Ty = VSUB(Tu, Tx);
						       T3h = VADD(Tu, Tx);
						       T1S = VZMULJ(Tf, T1R);
						       T1Q = VZMULJ(T1O, T1P);
						       T1I = VZMULJ(T7, T1H);
						       T3i = VADD(TB, TE);
						       TF = VSUB(TB, TE);
						       T1K = VZMULJ(T6, T1J);
						       T1A = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
						       T1E = LD(&(x[WS(rs, 19)]), ms, &(x[WS(rs, 1)]));
						  }
					     }
					     T35 = VADD(T1v, T1y);
					     T1z = VSUB(T1v, T1y);
					     T1T = VSUB(T1Q, T1S);
					     T36 = VADD(T1S, T1Q);
					     T3y = VSUB(T3h, T3i);
					     T3j = VADD(T3h, T3i);
					     T27 = VFMA(LDK(KP414213562), Ty, TF);
					     TG = VFNMS(LDK(KP414213562), TF, Ty);
					     T39 = VADD(T1I, T1K);
					     T1L = VSUB(T1I, T1K);
					     T1B = VZMULJ(T3, T1A);
					     T1F = VZMULJ(T1D, T1E);
					}
					{
					     V T11, T14, T1o, T1l, T1e, T1U, T1M, T1g, T16, T19;
					     {
						  V T10, T13, T1n, T1k;
						  T10 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
						  T13 = LD(&(x[WS(rs, 17)]), ms, &(x[WS(rs, 1)]));
						  T1n = LD(&(x[WS(rs, 25)]), ms, &(x[WS(rs, 1)]));
						  T1k = LD(&(x[WS(rs, 9)]), ms, &(x[WS(rs, 1)]));
						  {
						       V T1d, T1f, T1G, T38;
						       T1d = LD(&(x[WS(rs, 29)]), ms, &(x[WS(rs, 1)]));
						       T1f = LD(&(x[WS(rs, 13)]), ms, &(x[WS(rs, 1)]));
						       T1G = VSUB(T1B, T1F);
						       T38 = VADD(T1B, T1F);
						       T37 = VADD(T35, T36);
						       T3F = VSUB(T35, T36);
						       T11 = VZMULJ(T2, T10);
						       T14 = VZMULJ(T12, T13);
						       T1o = VZMULJ(T1m, T1n);
						       T1l = VZMULJ(T5, T1k);
						       T1e = VZMULJ(T1c, T1d);
						       T3G = VSUB(T39, T38);
						       T3a = VADD(T38, T39);
						       T1U = VSUB(T1L, T1G);
						       T1M = VADD(T1G, T1L);
						       T1g = VZMULJ(TK, T1f);
						  }
						  T16 = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
						  T19 = LD(&(x[WS(rs, 21)]), ms, &(x[WS(rs, 1)]));
					     }
					     T2Y = VADD(T11, T14);
					     T15 = VSUB(T11, T14);
					     T1p = VSUB(T1l, T1o);
					     T2Z = VADD(T1l, T1o);
					     T2w = VFNMS(LDK(KP707106781), T1U, T1T);
					     T1V = VFMA(LDK(KP707106781), T1U, T1T);
					     T2v = VFNMS(LDK(KP707106781), T1M, T1z);
					     T1N = VFMA(LDK(KP707106781), T1M, T1z);
					     T32 = VADD(T1e, T1g);
					     T1h = VSUB(T1e, T1g);
					     T17 = VZMULJ(TP, T16);
					     T1a = VZMULJ(T18, T19);
					}
				   }
			      }
			 }
			 {
			      V T2X, T3k, T3b, T3t, T1b, T31, T30, T3C, T3r, T3v, T3p, T3q;
			      T2X = VSUB(T2T, T2W);
			      T3p = VADD(T2T, T2W);
			      T3q = VADD(T3j, T3g);
			      T3k = VSUB(T3g, T3j);
			      T3b = VSUB(T37, T3a);
			      T3t = VADD(T37, T3a);
			      T1b = VSUB(T17, T1a);
			      T31 = VADD(T17, T1a);
			      T30 = VADD(T2Y, T2Z);
			      T3C = VSUB(T2Y, T2Z);
			      T3r = VADD(T3p, T3q);
			      T3v = VSUB(T3p, T3q);
			      {
				   V T3N, T3B, T3T, T3M, T3W, T3O, T2t, T1r, T2s, T1j, T3I, T3X, T3c, T3l, T3u;
				   V T3w;
				   {
					V T3L, T3A, T33, T3D, T1i, T1q;
					T3L = VSUB(T3z, T3y);
					T3A = VADD(T3y, T3z);
					T33 = VADD(T31, T32);
					T3D = VSUB(T31, T32);
					T1i = VADD(T1b, T1h);
					T1q = VSUB(T1b, T1h);
					{
					     V T3H, T3E, T34, T3s;
					     T3N = VFMA(LDK(KP414213562), T3F, T3G);
					     T3H = VFNMS(LDK(KP414213562), T3G, T3F);
					     T3B = VFMA(LDK(KP707106781), T3A, T3x);
					     T3T = VFNMS(LDK(KP707106781), T3A, T3x);
					     T3M = VFMA(LDK(KP707106781), T3L, T3K);
					     T3W = VFNMS(LDK(KP707106781), T3L, T3K);
					     T3O = VFMA(LDK(KP414213562), T3C, T3D);
					     T3E = VFNMS(LDK(KP414213562), T3D, T3C);
					     T34 = VSUB(T30, T33);
					     T3s = VADD(T30, T33);
					     T2t = VFNMS(LDK(KP707106781), T1q, T1p);
					     T1r = VFMA(LDK(KP707106781), T1q, T1p);
					     T2s = VFNMS(LDK(KP707106781), T1i, T15);
					     T1j = VFMA(LDK(KP707106781), T1i, T15);
					     T3I = VADD(T3E, T3H);
					     T3X = VSUB(T3H, T3E);
					     T3c = VADD(T34, T3b);
					     T3l = VSUB(T3b, T34);
					     T3u = VADD(T3s, T3t);
					     T3w = VSUB(T3t, T3s);
					}
				   }
				   {
					V T2p, Ts, TY, T1s, T2b, T2c, T1W, T26, T29, T2q, T3U, T3P, T2J, T2K;
					T2p = VFNMS(LDK(KP707106781), Tr, Tb);
					Ts = VFMA(LDK(KP707106781), Tr, Tb);
					T3U = VADD(T3O, T3N);
					T3P = VSUB(T3N, T3O);
					{
					     V T3Y, T40, T3R, T3J;
					     T3Y = VFMA(LDK(KP923879532), T3X, T3W);
					     T40 = VFNMS(LDK(KP923879532), T3X, T3W);
					     T3R = VFMA(LDK(KP923879532), T3I, T3B);
					     T3J = VFNMS(LDK(KP923879532), T3I, T3B);
					     {
						  V T3o, T3m, T3n, T3d;
						  T3o = VFMA(LDK(KP707106781), T3l, T3k);
						  T3m = VFNMS(LDK(KP707106781), T3l, T3k);
						  T3n = VFMA(LDK(KP707106781), T3c, T2X);
						  T3d = VFNMS(LDK(KP707106781), T3c, T2X);
						  ST(&(x[WS(rs, 24)]), VFNMSI(T3w, T3v), ms, &(x[0]));
						  ST(&(x[WS(rs, 8)]), VFMAI(T3w, T3v), ms, &(x[0]));
						  ST(&(x[0]), VADD(T3r, T3u), ms, &(x[0]));
						  ST(&(x[WS(rs, 16)]), VSUB(T3r, T3u), ms, &(x[0]));
						  {
						       V T3V, T3Z, T3S, T3Q;
						       T3V = VFNMS(LDK(KP923879532), T3U, T3T);
						       T3Z = VFMA(LDK(KP923879532), T3U, T3T);
						       T3S = VFMA(LDK(KP923879532), T3P, T3M);
						       T3Q = VFNMS(LDK(KP923879532), T3P, T3M);
						       ST(&(x[WS(rs, 4)]), VFMAI(T3o, T3n), ms, &(x[0]));
						       ST(&(x[WS(rs, 28)]), VFNMSI(T3o, T3n), ms, &(x[0]));
						       ST(&(x[WS(rs, 20)]), VFMAI(T3m, T3d), ms, &(x[0]));
						       ST(&(x[WS(rs, 12)]), VFNMSI(T3m, T3d), ms, &(x[0]));
						       ST(&(x[WS(rs, 22)]), VFNMSI(T3Y, T3V), ms, &(x[0]));
						       ST(&(x[WS(rs, 10)]), VFMAI(T3Y, T3V), ms, &(x[0]));
						       ST(&(x[WS(rs, 26)]), VFMAI(T40, T3Z), ms, &(x[0]));
						       ST(&(x[WS(rs, 6)]), VFNMSI(T40, T3Z), ms, &(x[0]));
						       ST(&(x[WS(rs, 2)]), VFMAI(T3S, T3R), ms, &(x[0]));
						       ST(&(x[WS(rs, 30)]), VFNMSI(T3S, T3R), ms, &(x[0]));
						       ST(&(x[WS(rs, 18)]), VFMAI(T3Q, T3J), ms, &(x[0]));
						       ST(&(x[WS(rs, 14)]), VFNMSI(T3Q, T3J), ms, &(x[0]));
						       TY = VADD(TG, TX);
						       T2B = VSUB(TX, TG);
						  }
					     }
					}
					T1s = VFNMS(LDK(KP198912367), T1r, T1j);
					T2b = VFMA(LDK(KP198912367), T1j, T1r);
					T2c = VFMA(LDK(KP198912367), T1N, T1V);
					T1W = VFNMS(LDK(KP198912367), T1V, T1N);
					T2A = VFMA(LDK(KP707106781), T25, T24);
					T26 = VFNMS(LDK(KP707106781), T25, T24);
					T29 = VSUB(T27, T28);
					T2q = VADD(T27, T28);
					{
					     V T2j, T2n, T1Y, T2f, T2o, T2m, T2e, T2g;
					     {
						  V T2h, TZ, T2i, T2d, T2l, T1X, T2k, T2a, T2D, T2E;
						  T2h = VFNMS(LDK(KP923879532), TY, Ts);
						  TZ = VFMA(LDK(KP923879532), TY, Ts);
						  T2i = VADD(T2b, T2c);
						  T2d = VSUB(T2b, T2c);
						  T2l = VSUB(T1W, T1s);
						  T1X = VADD(T1s, T1W);
						  T2k = VFNMS(LDK(KP923879532), T29, T26);
						  T2a = VFMA(LDK(KP923879532), T29, T26);
						  T2u = VFMA(LDK(KP668178637), T2t, T2s);
						  T2D = VFNMS(LDK(KP668178637), T2s, T2t);
						  T2j = VFNMS(LDK(KP980785280), T2i, T2h);
						  T2n = VFMA(LDK(KP980785280), T2i, T2h);
						  T2E = VFNMS(LDK(KP668178637), T2v, T2w);
						  T2x = VFMA(LDK(KP668178637), T2w, T2v);
						  T1Y = VFNMS(LDK(KP980785280), T1X, TZ);
						  T2f = VFMA(LDK(KP980785280), T1X, TZ);
						  T2o = VFMA(LDK(KP980785280), T2l, T2k);
						  T2m = VFNMS(LDK(KP980785280), T2l, T2k);
						  T2e = VFNMS(LDK(KP980785280), T2d, T2a);
						  T2g = VFMA(LDK(KP980785280), T2d, T2a);
						  T2r = VFMA(LDK(KP923879532), T2q, T2p);
						  T2J = VFNMS(LDK(KP923879532), T2q, T2p);
						  T2K = VADD(T2D, T2E);
						  T2F = VSUB(T2D, T2E);
					     }
					     ST(&(x[WS(rs, 23)]), VFMAI(T2m, T2j), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 9)]), VFNMSI(T2m, T2j), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 25)]), VFNMSI(T2o, T2n), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 7)]), VFMAI(T2o, T2n), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 31)]), VFMAI(T2g, T2f), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 1)]), VFNMSI(T2g, T2f), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 15)]), VFMAI(T2e, T1Y), ms, &(x[WS(rs, 1)]));
					     ST(&(x[WS(rs, 17)]), VFNMSI(T2e, T1Y), ms, &(x[WS(rs, 1)]));
					}
					T2L = VFMA(LDK(KP831469612), T2K, T2J);
					T2P = VFNMS(LDK(KP831469612), T2K, T2J);
				   }
			      }
			 }
		    }
	       }
	       {
		    V T2y, T2N, T2C, T2M;
		    T2y = VADD(T2u, T2x);
		    T2N = VSUB(T2x, T2u);
		    T2C = VFMA(LDK(KP923879532), T2B, T2A);
		    T2M = VFNMS(LDK(KP923879532), T2B, T2A);
		    {
			 V T2z, T2H, T2Q, T2O, T2G, T2I;
			 T2z = VFNMS(LDK(KP831469612), T2y, T2r);
			 T2H = VFMA(LDK(KP831469612), T2y, T2r);
			 T2Q = VFNMS(LDK(KP831469612), T2N, T2M);
			 T2O = VFMA(LDK(KP831469612), T2N, T2M);
			 T2G = VFNMS(LDK(KP831469612), T2F, T2C);
			 T2I = VFMA(LDK(KP831469612), T2F, T2C);
			 ST(&(x[WS(rs, 21)]), VFNMSI(T2O, T2L), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 11)]), VFMAI(T2O, T2L), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 27)]), VFMAI(T2Q, T2P), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 5)]), VFNMSI(T2Q, T2P), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 3)]), VFMAI(T2I, T2H), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 29)]), VFNMSI(T2I, T2H), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 19)]), VFMAI(T2G, T2z), ms, &(x[WS(rs, 1)]));
			 ST(&(x[WS(rs, 13)]), VFNMSI(T2G, T2z), ms, &(x[WS(rs, 1)]));
		    }
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 3),
     VTW(0, 9),
     VTW(0, 27),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 32, XSIMD_STRING("t3fv_32"), twinstr, &GENUS, {146, 116, 98, 0}, 0, 0, 0 };

void XSIMD(codelet_t3fv_32) (planner *p) {
     X(kdft_dit_register) (p, t3fv_32, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c.native -simd -compact -variables 4 -pipeline-latency 8 -twiddle-log3 -precompute-twiddles -no-generate-bytw -n 32 -name t3fv_32 -include t3f.h */

/*
 * This function contains 244 FP additions, 158 FP multiplications,
 * (or, 228 additions, 142 multiplications, 16 fused multiply/add),
 * 90 stack variables, 7 constants, and 64 memory accesses
 */
#include "t3f.h"

static void t3fv_32(R *ri, R *ii, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DVK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DVK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DVK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  R *x;
	  x = ri;
	  for (m = mb, W = W + (mb * ((TWVL / VL) * 8)); m < me; m = m + VL, x = x + (VL * ms), W = W + (TWVL * 8), MAKE_VOLATILE_STRIDE(rs)) {
	       V T2, T5, T3, T4, Tc, T1C, TP, Tz, Tn, T6, TS, Tf, TK, T7, T8;
	       V Tv, T1w, T22, Tg, Tk, T1D, T1R, TC, T18, T12, T1t, TH, TL, TT, T1n;
	       V T1c;
	       T2 = LDW(&(W[0]));
	       T5 = LDW(&(W[TWVL * 4]));
	       T3 = LDW(&(W[TWVL * 2]));
	       T4 = VZMULJ(T2, T3);
	       Tc = VZMUL(T2, T3);
	       T1C = VZMULJ(T2, T5);
	       TP = VZMULJ(T3, T5);
	       Tz = VZMUL(T2, T5);
	       Tn = VZMUL(T3, T5);
	       T6 = VZMUL(T4, T5);
	       TS = VZMULJ(Tc, T5);
	       Tf = VZMULJ(T4, T5);
	       TK = VZMUL(Tc, T5);
	       T7 = LDW(&(W[TWVL * 6]));
	       T8 = VZMULJ(T6, T7);
	       Tv = VZMULJ(T5, T7);
	       T1w = VZMULJ(Tn, T7);
	       T22 = VZMULJ(T3, T7);
	       Tg = VZMULJ(Tf, T7);
	       Tk = VZMUL(T2, T7);
	       T1D = VZMULJ(T1C, T7);
	       T1R = VZMULJ(Tc, T7);
	       TC = VZMULJ(T2, T7);
	       T18 = VZMULJ(TP, T7);
	       T12 = VZMULJ(Tz, T7);
	       T1t = VZMUL(Tc, T7);
	       TH = VZMUL(T3, T7);
	       TL = VZMULJ(TK, T7);
	       TT = VZMULJ(TS, T7);
	       T1n = VZMULJ(T4, T7);
	       T1c = VZMUL(T4, T7);
	       {
		    V Tb, T25, T2T, T3x, Tr, T1Z, T2W, T3K, TX, T27, T3g, T3z, TG, T28, T3j;
		    V T3y, T1N, T2v, T3a, T3G, T1V, T2w, T37, T3F, T1j, T2s, T33, T3D, T1r, T2t;
		    V T30, T3C;
		    {
			 V T1, T24, Ta, T21, T23, T9, T20, T2R, T2S;
			 T1 = LD(&(x[0]), ms, &(x[0]));
			 T23 = LD(&(x[WS(rs, 24)]), ms, &(x[0]));
			 T24 = VZMULJ(T22, T23);
			 T9 = LD(&(x[WS(rs, 16)]), ms, &(x[0]));
			 Ta = VZMULJ(T8, T9);
			 T20 = LD(&(x[WS(rs, 8)]), ms, &(x[0]));
			 T21 = VZMULJ(T1C, T20);
			 Tb = VSUB(T1, Ta);
			 T25 = VSUB(T21, T24);
			 T2R = VADD(T1, Ta);
			 T2S = VADD(T21, T24);
			 T2T = VADD(T2R, T2S);
			 T3x = VSUB(T2R, T2S);
		    }
		    {
			 V Te, Tp, Ti, Tm;
			 {
			      V Td, To, Th, Tl;
			      Td = LD(&(x[WS(rs, 4)]), ms, &(x[0]));
			      Te = VZMULJ(Tc, Td);
			      To = LD(&(x[WS(rs, 12)]), ms, &(x[0]));
			      Tp = VZMULJ(Tn, To);
			      Th = LD(&(x[WS(rs, 20)]), ms, &(x[0]));
			      Ti = VZMULJ(Tg, Th);
			      Tl = LD(&(x[WS(rs, 28)]), ms, &(x[0]));
			      Tm = VZMULJ(Tk, Tl);
			 }
			 {
			      V Tj, Tq, T2U, T2V;
			      Tj = VSUB(Te, Ti);
			      Tq = VSUB(Tm, Tp);
			      Tr = VMUL(LDK(KP707106781), VADD(Tj, Tq));
			      T1Z = VMUL(LDK(KP707106781), VSUB(Tq, Tj));
			      T2U = VADD(Te, Ti);
			      T2V = VADD(Tm, Tp);
			      T2W = VADD(T2U, T2V);
			      T3K = VSUB(T2V, T2U);
			 }
		    }
		    {
			 V TJ, TV, TN, TR;
			 {
			      V TI, TU, TM, TQ;
			      TI = LD(&(x[WS(rs, 30)]), ms, &(x[0]));
			      TJ = VZMULJ(TH, TI);
			      TU = LD(&(x[WS(rs, 22)]), ms, &(x[0]));
			      TV = VZMULJ(TT, TU);
			      TM = LD(&(x[WS(rs, 14)]), ms, &(x[0]));
			      TN = VZMULJ(TL, TM);
			      TQ = LD(&(x[WS(rs, 6)]), ms, &(x[0]));
			      TR = VZMULJ(TP, TQ);
			 }
			 {
			      V TO, TW, T3e, T3f;
			      TO = VSUB(TJ, TN);
			      TW = VSUB(TR, TV);
			      TX = VFMA(LDK(KP923879532), TO, VMUL(LDK(KP382683432), TW));
			      T27 = VFNMS(LDK(KP923879532), TW, VMUL(LDK(KP382683432), TO));
			      T3e = VADD(TJ, TN);
			      T3f = VADD(TR, TV);
			      T3g = VADD(T3e, T3f);
			      T3z = VSUB(T3e, T3f);
			 }
		    }
		    {
			 V Tu, TE, Tx, TB;
			 {
			      V Tt, TD, Tw, TA;
			      Tt = LD(&(x[WS(rs, 2)]), ms, &(x[0]));
			      Tu = VZMULJ(T4, Tt);
			      TD = LD(&(x[WS(rs, 26)]), ms, &(x[0]));
			      TE = VZMULJ(TC, TD);
			      Tw = LD(&(x[WS(rs, 18)]), ms, &(x[0]));
			      Tx = VZMULJ(Tv, Tw);
			      TA = LD(&(x[WS(rs, 10)]), ms, &(x[0]));
			      TB = VZMULJ(Tz, TA);
			 }
			 {
			      V Ty, TF, T3h, T3i;
			      Ty = VSUB(Tu, Tx);
			      TF = VSUB(TB, TE);
			      TG = VFNMS(LDK(KP382683432), TF, VMUL(LDK(KP923879532), Ty));
			      T28 = VFMA(LDK(KP382683432), Ty, VMUL(LDK(KP923879532), TF));
			      T3h = VADD(Tu, Tx);
			      T3i = VADD(TB, TE);
			      T3j = VADD(T3h, T3i);
			      T3y = VSUB(T3h, T3i);
			 }
		    }
		    {
			 V T1v, T1y, T1T, T1Q, T1I, T1K, T1L, T1B, T1F, T1G;
			 {
			      V T1u, T1x, T1S, T1P;
			      T1u = LD(&(x[WS(rs, 31)]), ms, &(x[WS(rs, 1)]));
			      T1v = VZMULJ(T1t, T1u);
			      T1x = LD(&(x[WS(rs, 15)]), ms, &(x[WS(rs, 1)]));
			      T1y = VZMULJ(T1w, T1x);
			      T1S = LD(&(x[WS(rs, 23)]), ms, &(x[WS(rs, 1)]));
			      T1T = VZMULJ(T1R, T1S);
			      T1P = LD(&(x[WS(rs, 7)]), ms, &(x[WS(rs, 1)]));
			      T1Q = VZMULJ(Tf, T1P);
			      {
				   V T1H, T1J, T1A, T1E;
				   T1H = LD(&(x[WS(rs, 27)]), ms, &(x[WS(rs, 1)]));
				   T1I = VZMULJ(T7, T1H);
				   T1J = LD(&(x[WS(rs, 11)]), ms, &(x[WS(rs, 1)]));
				   T1K = VZMULJ(T6, T1J);
				   T1L = VSUB(T1I, T1K);
				   T1A = LD(&(x[WS(rs, 3)]), ms, &(x[WS(rs, 1)]));
				   T1B = VZMULJ(T3, T1A);
				   T1E = LD(&(x[WS(rs, 19)]), ms, &(x[WS(rs, 1)]));
				   T1F = VZMULJ(T1D, T1E);
				   T1G = VSUB(T1B, T1F);
			      }
			 }
			 {
			      V T1z, T1M, T38, T39;
			      T1z = VSUB(T1v, T1y);
			      T1M = VMUL(LDK(KP707106781), VADD(T1G, T1L));
			      T1N = VADD(T1z, T1M);
			      T2v = VSUB(T1z, T1M);
			      T38 = VADD(T1B, T1F);
			      T39 = VADD(T1I, T1K);
			      T3a = VADD(T38, T39);
			      T3G = VSUB(T39, T38);
			 }
			 {
			      V T1O, T1U, T35, T36;
			      T1O = VMUL(LDK(KP707106781), VSUB(T1L, T1G));
			      T1U = VSUB(T1Q, T1T);
			      T1V = VSUB(T1O, T1U);
			      T2w = VADD(T1U, T1O);
			      T35 = VADD(T1v, T1y);
			      T36 = VADD(T1Q, T1T);
			      T37 = VADD(T35, T36);
			      T3F = VSUB(T35, T36);
			 }
		    }
		    {
			 V T11, T14, T1p, T1m, T1e, T1g, T1h, T17, T1a, T1b;
			 {
			      V T10, T13, T1o, T1l;
			      T10 = LD(&(x[WS(rs, 1)]), ms, &(x[WS(rs, 1)]));
			      T11 = VZMULJ(T2, T10);
			      T13 = LD(&(x[WS(rs, 17)]), ms, &(x[WS(rs, 1)]));
			      T14 = VZMULJ(T12, T13);
			      T1o = LD(&(x[WS(rs, 25)]), ms, &(x[WS(rs, 1)]));
			      T1p = VZMULJ(T1n, T1o);
			      T1l = LD(&(x[WS(rs, 9)]), ms, &(x[WS(rs, 1)]));
			      T1m = VZMULJ(T5, T1l);
			      {
				   V T1d, T1f, T16, T19;
				   T1d = LD(&(x[WS(rs, 29)]), ms, &(x[WS(rs, 1)]));
				   T1e = VZMULJ(T1c, T1d);
				   T1f = LD(&(x[WS(rs, 13)]), ms, &(x[WS(rs, 1)]));
				   T1g = VZMULJ(TK, T1f);
				   T1h = VSUB(T1e, T1g);
				   T16 = LD(&(x[WS(rs, 5)]), ms, &(x[WS(rs, 1)]));
				   T17 = VZMULJ(TS, T16);
				   T19 = LD(&(x[WS(rs, 21)]), ms, &(x[WS(rs, 1)]));
				   T1a = VZMULJ(T18, T19);
				   T1b = VSUB(T17, T1a);
			      }
			 }
			 {
			      V T15, T1i, T31, T32;
			      T15 = VSUB(T11, T14);
			      T1i = VMUL(LDK(KP707106781), VADD(T1b, T1h));
			      T1j = VADD(T15, T1i);
			      T2s = VSUB(T15, T1i);
			      T31 = VADD(T17, T1a);
			      T32 = VADD(T1e, T1g);
			      T33 = VADD(T31, T32);
			      T3D = VSUB(T32, T31);
			 }
			 {
			      V T1k, T1q, T2Y, T2Z;
			      T1k = VMUL(LDK(KP707106781), VSUB(T1h, T1b));
			      T1q = VSUB(T1m, T1p);
			      T1r = VSUB(T1k, T1q);
			      T2t = VADD(T1q, T1k);
			      T2Y = VADD(T11, T14);
			      T2Z = VADD(T1m, T1p);
			      T30 = VADD(T2Y, T2Z);
			      T3C = VSUB(T2Y, T2Z);
			 }
		    }
		    {
			 V T3r, T3v, T3u, T3w;
			 {
			      V T3p, T3q, T3s, T3t;
			      T3p = VADD(T2T, T2W);
			      T3q = VADD(T3j, T3g);
			      T3r = VADD(T3p, T3q);
			      T3v = VSUB(T3p, T3q);
			      T3s = VADD(T30, T33);
			      T3t = VADD(T37, T3a);
			      T3u = VADD(T3s, T3t);
			      T3w = VBYI(VSUB(T3t, T3s));
			 }
			 ST(&(x[WS(rs, 16)]), VSUB(T3r, T3u), ms, &(x[0]));
			 ST(&(x[WS(rs, 8)]), VADD(T3v, T3w), ms, &(x[0]));
			 ST(&(x[0]), VADD(T3r, T3u), ms, &(x[0]));
			 ST(&(x[WS(rs, 24)]), VSUB(T3v, T3w), ms, &(x[0]));
		    }
		    {
			 V T2X, T3k, T3c, T3l, T34, T3b;
			 T2X = VSUB(T2T, T2W);
			 T3k = VSUB(T3g, T3j);
			 T34 = VSUB(T30, T33);
			 T3b = VSUB(T37, T3a);
			 T3c = VMUL(LDK(KP707106781), VADD(T34, T3b));
			 T3l = VMUL(LDK(KP707106781), VSUB(T3b, T34));
			 {
			      V T3d, T3m, T3n, T3o;
			      T3d = VADD(T2X, T3c);
			      T3m = VBYI(VADD(T3k, T3l));
			      ST(&(x[WS(rs, 28)]), VSUB(T3d, T3m), ms, &(x[0]));
			      ST(&(x[WS(rs, 4)]), VADD(T3d, T3m), ms, &(x[0]));
			      T3n = VSUB(T2X, T3c);
			      T3o = VBYI(VSUB(T3l, T3k));
			      ST(&(x[WS(rs, 20)]), VSUB(T3n, T3o), ms, &(x[0]));
			      ST(&(x[WS(rs, 12)]), VADD(T3n, T3o), ms, &(x[0]));
			 }
		    }
		    {
			 V T3B, T3W, T3M, T3U, T3I, T3T, T3P, T3X, T3A, T3L;
			 T3A = VMUL(LDK(KP707106781), VADD(T3y, T3z));
			 T3B = VADD(T3x, T3A);
			 T3W = VSUB(T3x, T3A);
			 T3L = VMUL(LDK(KP707106781), VSUB(T3z, T3y));
			 T3M = VADD(T3K, T3L);
			 T3U = VSUB(T3L, T3K);
			 {
			      V T3E, T3H, T3N, T3O;
			      T3E = VFMA(LDK(KP923879532), T3C, VMUL(LDK(KP382683432), T3D));
			      T3H = VFNMS(LDK(KP382683432), T3G, VMUL(LDK(KP923879532), T3F));
			      T3I = VADD(T3E, T3H);
			      T3T = VSUB(T3H, T3E);
			      T3N = VFNMS(LDK(KP382683432), T3C, VMUL(LDK(KP923879532), T3D));
			      T3O = VFMA(LDK(KP382683432), T3F, VMUL(LDK(KP923879532), T3G));
			      T3P = VADD(T3N, T3O);
			      T3X = VSUB(T3O, T3N);
			 }
			 {
			      V T3J, T3Q, T3Z, T40;
			      T3J = VADD(T3B, T3I);
			      T3Q = VBYI(VADD(T3M, T3P));
			      ST(&(x[WS(rs, 30)]), VSUB(T3J, T3Q), ms, &(x[0]));
			      ST(&(x[WS(rs, 2)]), VADD(T3J, T3Q), ms, &(x[0]));
			      T3Z = VBYI(VADD(T3U, T3T));
			      T40 = VADD(T3W, T3X);
			      ST(&(x[WS(rs, 6)]), VADD(T3Z, T40), ms, &(x[0]));
			      ST(&(x[WS(rs, 26)]), VSUB(T40, T3Z), ms, &(x[0]));
			 }
			 {
			      V T3R, T3S, T3V, T3Y;
			      T3R = VSUB(T3B, T3I);
			      T3S = VBYI(VSUB(T3P, T3M));
			      ST(&(x[WS(rs, 18)]), VSUB(T3R, T3S), ms, &(x[0]));
			      ST(&(x[WS(rs, 14)]), VADD(T3R, T3S), ms, &(x[0]));
			      T3V = VBYI(VSUB(T3T, T3U));
			      T3Y = VSUB(T3W, T3X);
			      ST(&(x[WS(rs, 10)]), VADD(T3V, T3Y), ms, &(x[0]));
			      ST(&(x[WS(rs, 22)]), VSUB(T3Y, T3V), ms, &(x[0]));
			 }
		    }
		    {
			 V TZ, T2k, T2d, T2l, T1X, T2h, T2a, T2i;
			 {
			      V Ts, TY, T2b, T2c;
			      Ts = VADD(Tb, Tr);
			      TY = VADD(TG, TX);
			      TZ = VADD(Ts, TY);
			      T2k = VSUB(Ts, TY);
			      T2b = VFNMS(LDK(KP195090322), T1j, VMUL(LDK(KP980785280), T1r));
			      T2c = VFMA(LDK(KP195090322), T1N, VMUL(LDK(KP980785280), T1V));
			      T2d = VADD(T2b, T2c);
			      T2l = VSUB(T2c, T2b);
			 }
			 {
			      V T1s, T1W, T26, T29;
			      T1s = VFMA(LDK(KP980785280), T1j, VMUL(LDK(KP195090322), T1r));
			      T1W = VFNMS(LDK(KP195090322), T1V, VMUL(LDK(KP980785280), T1N));
			      T1X = VADD(T1s, T1W);
			      T2h = VSUB(T1W, T1s);
			      T26 = VSUB(T1Z, T25);
			      T29 = VSUB(T27, T28);
			      T2a = VADD(T26, T29);
			      T2i = VSUB(T29, T26);
			 }
			 {
			      V T1Y, T2e, T2n, T2o;
			      T1Y = VADD(TZ, T1X);
			      T2e = VBYI(VADD(T2a, T2d));
			      ST(&(x[WS(rs, 31)]), VSUB(T1Y, T2e), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 1)]), VADD(T1Y, T2e), ms, &(x[WS(rs, 1)]));
			      T2n = VBYI(VADD(T2i, T2h));
			      T2o = VADD(T2k, T2l);
			      ST(&(x[WS(rs, 7)]), VADD(T2n, T2o), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 25)]), VSUB(T2o, T2n), ms, &(x[WS(rs, 1)]));
			 }
			 {
			      V T2f, T2g, T2j, T2m;
			      T2f = VSUB(TZ, T1X);
			      T2g = VBYI(VSUB(T2d, T2a));
			      ST(&(x[WS(rs, 17)]), VSUB(T2f, T2g), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 15)]), VADD(T2f, T2g), ms, &(x[WS(rs, 1)]));
			      T2j = VBYI(VSUB(T2h, T2i));
			      T2m = VSUB(T2k, T2l);
			      ST(&(x[WS(rs, 9)]), VADD(T2j, T2m), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 23)]), VSUB(T2m, T2j), ms, &(x[WS(rs, 1)]));
			 }
		    }
		    {
			 V T2r, T2M, T2F, T2N, T2y, T2J, T2C, T2K;
			 {
			      V T2p, T2q, T2D, T2E;
			      T2p = VSUB(Tb, Tr);
			      T2q = VADD(T28, T27);
			      T2r = VADD(T2p, T2q);
			      T2M = VSUB(T2p, T2q);
			      T2D = VFNMS(LDK(KP555570233), T2s, VMUL(LDK(KP831469612), T2t));
			      T2E = VFMA(LDK(KP555570233), T2v, VMUL(LDK(KP831469612), T2w));
			      T2F = VADD(T2D, T2E);
			      T2N = VSUB(T2E, T2D);
			 }
			 {
			      V T2u, T2x, T2A, T2B;
			      T2u = VFMA(LDK(KP831469612), T2s, VMUL(LDK(KP555570233), T2t));
			      T2x = VFNMS(LDK(KP555570233), T2w, VMUL(LDK(KP831469612), T2v));
			      T2y = VADD(T2u, T2x);
			      T2J = VSUB(T2x, T2u);
			      T2A = VADD(T25, T1Z);
			      T2B = VSUB(TX, TG);
			      T2C = VADD(T2A, T2B);
			      T2K = VSUB(T2B, T2A);
			 }
			 {
			      V T2z, T2G, T2P, T2Q;
			      T2z = VADD(T2r, T2y);
			      T2G = VBYI(VADD(T2C, T2F));
			      ST(&(x[WS(rs, 29)]), VSUB(T2z, T2G), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 3)]), VADD(T2z, T2G), ms, &(x[WS(rs, 1)]));
			      T2P = VBYI(VADD(T2K, T2J));
			      T2Q = VADD(T2M, T2N);
			      ST(&(x[WS(rs, 5)]), VADD(T2P, T2Q), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 27)]), VSUB(T2Q, T2P), ms, &(x[WS(rs, 1)]));
			 }
			 {
			      V T2H, T2I, T2L, T2O;
			      T2H = VSUB(T2r, T2y);
			      T2I = VBYI(VSUB(T2F, T2C));
			      ST(&(x[WS(rs, 19)]), VSUB(T2H, T2I), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 13)]), VADD(T2H, T2I), ms, &(x[WS(rs, 1)]));
			      T2L = VBYI(VSUB(T2J, T2K));
			      T2O = VSUB(T2M, T2N);
			      ST(&(x[WS(rs, 11)]), VADD(T2L, T2O), ms, &(x[WS(rs, 1)]));
			      ST(&(x[WS(rs, 21)]), VSUB(T2O, T2L), ms, &(x[WS(rs, 1)]));
			 }
		    }
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(0, 1),
     VTW(0, 3),
     VTW(0, 9),
     VTW(0, 27),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 32, XSIMD_STRING("t3fv_32"), twinstr, &GENUS, {228, 142, 16, 0}, 0, 0, 0 };

void XSIMD(codelet_t3fv_32) (planner *p) {
     X(kdft_dit_register) (p, t3fv_32, &desc);
}
#endif				/* HAVE_FMA */
