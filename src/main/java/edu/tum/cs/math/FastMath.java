/* Adapted from the jafama library written by Jeff Hain (https://github.com/jeffhain/jafama).
   We would like to simply import the library as is, but it contains some optimizations than can only be enabled via
   Java system properties (jafama.fastlog). Since those can only be set at JVM startup, we found it easier and more
   reliable to import the relevant parts of the code. */
/*
 * Copyright 2012-2015 Jeff Hain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * =============================================================================
 * Notice of fdlibm package this program is partially derived from:
 *
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * =============================================================================
 */
package edu.tum.cs.math;

public class FastMath {

	static final double DOUBLE_MIN_NORMAL = Double.longBitsToDouble(0x0010000000000000L); // 2.2250738585072014E-308

    // Not storing float/double mantissa size in constants,
    // for 23 and 52 are shorter to read and more
    // bitwise-explicit than some constant's name.

    static final int MIN_DOUBLE_EXPONENT = -1074;
    static final int MAX_DOUBLE_EXPONENT = 1023;

    static final double EXP_OVERFLOW_LIMIT = Double.longBitsToDouble(0x40862E42FEFA39EFL); // 7.09782712893383973096e+02
    static final double EXP_UNDERFLOW_LIMIT = Double.longBitsToDouble(0xC0874910D52D3051L); // -7.45133219101941108420e+02
    static final int EXP_LO_DISTANCE_TO_ZERO_POT = 0;
    static final int EXP_LO_DISTANCE_TO_ZERO = (1<<EXP_LO_DISTANCE_TO_ZERO_POT);
    static final int EXP_LO_TAB_SIZE_POT = 11;
    static final int EXP_LO_TAB_SIZE = (1<<EXP_LO_TAB_SIZE_POT)+1;
    static final int EXP_LO_TAB_MID_INDEX = ((EXP_LO_TAB_SIZE-1)/2);
    static final int EXP_LO_INDEXING = EXP_LO_TAB_MID_INDEX/EXP_LO_DISTANCE_TO_ZERO;
    static final int EXP_LO_INDEXING_DIV_SHIFT = EXP_LO_TAB_SIZE_POT-1-EXP_LO_DISTANCE_TO_ZERO_POT;
    
    static final class MyTExp {
        static final double[] expHiTab = new double[1+(int)EXP_OVERFLOW_LIMIT-(int)EXP_UNDERFLOW_LIMIT];
        static final double[] expLoPosTab = new double[EXP_LO_TAB_SIZE];
        static final double[] expLoNegTab = new double[EXP_LO_TAB_SIZE];
        static {
            init();
        }
        private static strictfp void init() {
            for (int i=(int)EXP_UNDERFLOW_LIMIT;i<=(int)EXP_OVERFLOW_LIMIT;i++) {
                expHiTab[i-(int)EXP_UNDERFLOW_LIMIT] = StrictMath.exp(i);
            }
            for (int i=0;i<EXP_LO_TAB_SIZE;i++) {
                // x: in [-EXPM1_DISTANCE_TO_ZERO,EXPM1_DISTANCE_TO_ZERO].
                double x = -EXP_LO_DISTANCE_TO_ZERO + i/(double)EXP_LO_INDEXING;
                // exp(x)
                expLoPosTab[i] = StrictMath.exp(x);
                // 1-exp(-x), accurately computed
                expLoNegTab[i] = -StrictMath.expm1(-x);
            }
        }
    }

    static final int LOG_BITS = 12;
    static final int LOG_TAB_SIZE = (1<<LOG_BITS);
    
    static final class MyTLog {
        static final double[] logXLogTab = new double[LOG_TAB_SIZE];
        static final double[] logXTab = new double[LOG_TAB_SIZE];
        static final double[] logXInvTab = new double[LOG_TAB_SIZE];
        static {
            init();
        }
        private static strictfp void init() {
            for (int i=0;i<LOG_TAB_SIZE;i++) {
                // Exact to use inverse of tab size, since it is a power of two.
                double x = 1+i*(1.0/LOG_TAB_SIZE);
                logXLogTab[i] = StrictMath.log(x);
                logXTab[i] = x;
                logXInvTab[i] = 1/x;
            }
        }
    }

    static final boolean USE_TWO_POW_TAB = true;
    static final int TWO_POW_TAB_SIZE = USE_TWO_POW_TAB ? (MAX_DOUBLE_EXPONENT-MIN_DOUBLE_EXPONENT)+1 : 0;
    
    static final class MyTTwoPow {
        static final double[] twoPowTab = new double[TWO_POW_TAB_SIZE];
        static {
            init();
        }
        private static strictfp void init() {
            if (USE_TWO_POW_TAB) {
                for (int i=MIN_DOUBLE_EXPONENT;i<=MAX_DOUBLE_EXPONENT;i++) {
                    twoPowTab[i-MIN_DOUBLE_EXPONENT] = twoPow(i);
                }
            }
        }
    }

    /**
     * Returns the exact result, provided it's in double range,
     * i.e. if power is in [-1074,1023].
     * 
     * @param power An int power.
     * @return 2^power as a double, or +-Infinity in case of overflow.
     */
    public static double twoPow(int power) {
        if (power <= -MAX_DOUBLE_EXPONENT) { // Not normal.
            if (power >= MIN_DOUBLE_EXPONENT) { // Subnormal.
                return Double.longBitsToDouble(0x0008000000000000L>>(-(power+MAX_DOUBLE_EXPONENT)));
            } else { // Underflow.
                return 0.0;
            }
        } else if (power > MAX_DOUBLE_EXPONENT) { // Overflow.
            return Double.POSITIVE_INFINITY;
        } else { // Normal.
            return Double.longBitsToDouble(((long)(power+MAX_DOUBLE_EXPONENT))<<52);
        }
    }

    /**
     * @param power Must be in normal or subnormal values range.
     */
    static double twoPowNormalOrSubnormal(int power) {
        if (USE_TWO_POW_TAB) {
            return MyTTwoPow.twoPowTab[power-MIN_DOUBLE_EXPONENT];
        } else {
            if (power <= -MAX_DOUBLE_EXPONENT) { // Not normal.
                return Double.longBitsToDouble(0x0008000000000000L>>(-(power+MAX_DOUBLE_EXPONENT)));
            } else { // Normal.
                return Double.longBitsToDouble(((long)(power+MAX_DOUBLE_EXPONENT))<<52);
            }
        }
    }

    static final double TWO_POW_52 = twoPow(52);

    static final double LOG_2 = StrictMath.log(2.0);

    /**
     * @param value A double value.
     * @return Value logarithm (base e).
     */
    public static double log(double value) {
        if (value > 0.0) {
            if (value == Double.POSITIVE_INFINITY) {
                return Double.POSITIVE_INFINITY;
            }

            // For normal values not close to 1.0, we use the following formula:
            // log(value)
            // = log(2^exponent*1.mantissa)
            // = log(2^exponent) + log(1.mantissa)
            // = exponent * log(2) + log(1.mantissa)
            // = exponent * log(2) + log(1.mantissaApprox) + log(1.mantissa/1.mantissaApprox)
            // = exponent * log(2) + log(1.mantissaApprox) + log(1+epsilon)
            // = exponent * log(2) + log(1.mantissaApprox) + epsilon-epsilon^2/2+epsilon^3/3-epsilon^4/4+...
            // with:
            // 1.mantissaApprox <= 1.mantissa,
            // log(1.mantissaApprox) in table,
            // epsilon = (1.mantissa/1.mantissaApprox)-1
            //
            // To avoid bad relative error for small results,
            // values close to 1.0 are treated aside, with the formula:
            // log(x) = z*(2+z^2*((2.0/3)+z^2*((2.0/5))+z^2*((2.0/7))+...)))
            // with z=(x-1)/(x+1)

            double h;
            if (value > 0.95) {
                if (value < 1.14) {
                    double z = (value-1.0)/(value+1.0);
                    double z2 = z*z;
                    return z*(2+z2*((2.0/3)+z2*((2.0/5)+z2*((2.0/7)+z2*((2.0/9)+z2*((2.0/11)))))));
                }
                h = 0.0;
            } else if (value < DOUBLE_MIN_NORMAL) {
                // Ensuring value is normal.
                value *= TWO_POW_52;
                // log(x*2^52)
                // = log(x)-ln(2^52)
                // = log(x)-52*ln(2)
                h = -52*LOG_2;
            } else {
                h = 0.0;
            }

            int valueBitsHi = (int)(Double.doubleToRawLongBits(value)>>32);
            int valueExp = (valueBitsHi>>20)-MAX_DOUBLE_EXPONENT;
            // Getting the first LOG_BITS bits of the mantissa.
            int xIndex = ((valueBitsHi<<12)>>>(32-LOG_BITS));

            // 1.mantissa/1.mantissaApprox - 1
            double z = (value * twoPowNormalOrSubnormal(-valueExp)) * MyTLog.logXInvTab[xIndex] - 1;

            z *= (1-z*((1.0/2)-z*((1.0/3))));

            return h + valueExp * LOG_2 + (MyTLog.logXLogTab[xIndex] + z);

        } else if (value == 0.0) {
            return Double.NEGATIVE_INFINITY;
        } else { // value < 0.0, or value is NaN
            return Double.NaN;
        }
    }

    /**
     * @param value A double value.
     * @return e^value.
     */
    public static double exp(double value) {
        // exp(x) = exp([x])*exp(y)
        // with [x] the integer part of x, and y = x-[x]
        // ===>
        // We find an approximation of y, called z.
        // ===>
        // exp(x) = exp([x])*(exp(z)*exp(epsilon))
        // with epsilon = y - z
        // ===>
        // We have exp([x]) and exp(z) pre-computed in tables, we "just" have to compute exp(epsilon).
        //
        // We use the same indexing (cast to int) to compute x integer part and the
        // table index corresponding to z, to avoid two int casts.
        // Also, to optimize index multiplication and division, we use powers of two,
        // so that we can do it with bits shifts.

        if (value > EXP_OVERFLOW_LIMIT) {
            return Double.POSITIVE_INFINITY;
        } else if (!(value >= EXP_UNDERFLOW_LIMIT)) {
            return (value != value) ? Double.NaN : 0.0;
        }

        final int indexes = (int)(value*EXP_LO_INDEXING);

        final int valueInt;
        if (indexes >= 0) {
            valueInt = (indexes>>EXP_LO_INDEXING_DIV_SHIFT);
        } else {
            valueInt = -((-indexes)>>EXP_LO_INDEXING_DIV_SHIFT);
        }
        final double hiTerm = MyTExp.expHiTab[valueInt-(int)EXP_UNDERFLOW_LIMIT];

        final int zIndex = indexes - (valueInt<<EXP_LO_INDEXING_DIV_SHIFT);
        final double y = (value-valueInt);
        final double z = zIndex*(1.0/EXP_LO_INDEXING);
        final double eps = y-z;
        final double expZ = MyTExp.expLoPosTab[zIndex+EXP_LO_TAB_MID_INDEX];
        final double expEps = (1+eps*(1+eps*(1.0/2+eps*(1.0/6+eps*(1.0/24)))));
        final double loTerm = expZ * expEps;

        return hiTerm * loTerm;
    }

    /**
     * 1e-13ish accuracy or better on whole double range.
     * 
     * @param value A double value.
     * @param power A power.
     * @return value^power.
     */
    public static double pow(double value, double power) {
        if (power == 0.0) {
            return 1.0;
        } else if (power == 1.0) {
            return value;
        }
        if (value <= 0.0) {
            // powerInfo: 0 if not integer, 1 if even integer, -1 if odd integer
            int powerInfo;
            if (Math.abs(power) >= (TWO_POW_52*2)) {
                // The binary digit just before comma is outside mantissa,
                // thus it is always 0: power is an even integer.
                powerInfo = 1;
            } else {
                // If power's magnitude permits, we cast into int instead of into long,
                // as it is faster.
                if (Math.abs(power) <= (double)Integer.MAX_VALUE) {
                    int powerAsInt = (int)power;
                    if (power == (double)powerAsInt) {
                        powerInfo = ((powerAsInt & 1) == 0) ? 1 : -1;
                    } else { // power is not an integer (and not NaN, due to test against Integer.MAX_VALUE)
                        powerInfo = 0;
                    }
                } else {
                    long powerAsLong = (long)power;
                    if (power == (double)powerAsLong) {
                        powerInfo = ((powerAsLong & 1) == 0) ? 1 : -1;
                    } else { // power is not an integer, or is NaN
                        if (power != power) {
                            return Double.NaN;
                        }
                        powerInfo = 0;
                    }
                }
            }

            if (value == 0.0) {
                if (power < 0.0) {
                    return (powerInfo < 0) ? 1/value : Double.POSITIVE_INFINITY;
                } else { // power > 0.0 (0 and NaN cases already treated)
                    return (powerInfo < 0) ? value : 0.0;
                }
            } else { // value < 0.0
                if (value == Double.NEGATIVE_INFINITY) {
                    if (powerInfo < 0) { // power odd integer
                        return (power < 0.0) ? -0.0 : Double.NEGATIVE_INFINITY;
                    } else { // power even integer, or not an integer
                        return (power < 0.0) ? 0.0 : Double.POSITIVE_INFINITY;
                    }
                } else {
                    return (powerInfo == 0) ? Double.NaN : powerInfo * exp(power*log(-value));
                }
            }
        } else { // value > 0.0, or value is NaN
            return exp(power*log(value));
        }
    }

    public static double log1p(double value) {
        if (value > -1.0) {
            if (value == Double.POSITIVE_INFINITY) {
                return Double.POSITIVE_INFINITY;
            }

            // ln'(x) = 1/x
            // so
            // log(x+epsilon) ~= log(x) + epsilon/x
            // 
            // Let u be 1+value rounded:
            // 1+value = u+epsilon
            //
            // log(1+value)
            // = log(u+epsilon)
            // ~= log(u) + epsilon/value
            // We compute log(u) as done in log(double), and then add the corrective term.

            double valuePlusOne = 1.0+value;
            if (valuePlusOne == 1.0) {
                return value;
            } else if (Math.abs(value) < 0.15) {
                double z = value/(value+2.0);
                double z2 = z*z;
                return z*(2+z2*((2.0/3)+z2*((2.0/5)+z2*((2.0/7)+z2*((2.0/9)+z2*((2.0/11)))))));
            }

            int valuePlusOneBitsHi = (int)(Double.doubleToRawLongBits(valuePlusOne)>>32) & 0x7FFFFFFF;
            int valuePlusOneExp = (valuePlusOneBitsHi>>20)-MAX_DOUBLE_EXPONENT;
            // Getting the first LOG_BITS bits of the mantissa.
            int xIndex = ((valuePlusOneBitsHi<<12)>>>(32-LOG_BITS));

            // 1.mantissa/1.mantissaApprox - 1
            double z = (valuePlusOne * twoPowNormalOrSubnormal(-valuePlusOneExp)) * MyTLog.logXInvTab[xIndex] - 1;

            z *= (1-z*((1.0/2)-z*(1.0/3)));

            // Adding epsilon/valuePlusOne to z,
            // with
            // epsilon = value - (valuePlusOne-1)
            // (valuePlusOne + epsilon ~= 1+value (not rounded))

            return valuePlusOneExp * LOG_2 + MyTLog.logXLogTab[xIndex] + (z + (value - (valuePlusOne-1))/valuePlusOne);
        } else if (value == -1.0) {
            return Double.NEGATIVE_INFINITY;
        } else { // value < -1.0, or value is NaN
            return Double.NaN;
        }
    }

}
