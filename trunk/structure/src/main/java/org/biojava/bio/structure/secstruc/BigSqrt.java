/*
 *                    PDB web development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 *
 * Created on Aug 16, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.bio.structure.secstruc;

import java.math.BigDecimal;
import java.math.BigInteger;

/** calculate a more precise SQRT.
 *
 * Modified from:
 * http://www.merriampark.com/bigsqrt.htm
 *
 * @author Andreas Prlic
 *
 */
public class BigSqrt {

   private static BigDecimal ZERO = new BigDecimal ("0");
   private static BigDecimal ONE = new BigDecimal ("1");
   private static BigDecimal TWO = new BigDecimal ("2");

   public static final int DEFAULT_MAX_ITERATIONS = 50;

   /** we take 3 aftercomma
    *
    */
   public static final int DEFAULT_SCALE = 3;

   private BigDecimal error;
   private int iterations;
   private boolean traceFlag;
   private int scale = DEFAULT_SCALE;
   private int maxIterations = DEFAULT_MAX_ITERATIONS;

   //---------------------------------------
   // The error is the original number minus
   // (sqrt * sqrt). If the original number
   // was a perfect square, the error is 0.
   //---------------------------------------

   public BigDecimal getError () {
     return error;
   }

   //-------------------------------------------------------------
   // Number of iterations performed when square root was computed
   //-------------------------------------------------------------

   public int getIterations () {
     return iterations;
   }

   //-----------
   // Trace flag
   //-----------

   public boolean getTraceFlag () {
     return traceFlag;
   }

   public void setTraceFlag (boolean flag) {
     traceFlag = flag;
   }

   //------
   // Scale
   //------

   public int getScale () {
     return scale;
   }

   public void setScale (int scale) {
     this.scale = scale;
   }

   //-------------------
   // Maximum iterations
   //-------------------

   public int getMaxIterations () {
     return maxIterations;
   }

   public void setMaxIterations (int maxIterations) {
     this.maxIterations = maxIterations;
   }

   //--------------------------
   // Get initial approximation
   //--------------------------

   private static BigDecimal getInitialApproximation (BigDecimal n) {
     BigInteger integerPart = n.toBigInteger ();
     int length = integerPart.toString ().length ();
     if ((length % 2) == 0) {
       length--;
     }
     length /= 2;
     BigDecimal guess = ONE.movePointRight (length);
     return guess;
   }

   /** Get square root
    *
    * @param n
    * @return a BigDecimal
    */


   public BigDecimal sqrt (BigInteger n) {
     return sqrt (new BigDecimal (n));
   }

   /** Get square root
   *
   * @param n
   * @return a BigDecimal
   */
   public BigDecimal sqrt(BigDecimal n) {

     // Make sure n is a positive number

     if (n.compareTo (ZERO) <= 0) {
       throw new IllegalArgumentException ();
     }

     BigDecimal initialGuess = getInitialApproximation (n);
     //System.out.println("Initial guess " + initialGuess.toString ());
     BigDecimal lastGuess = ZERO;
     BigDecimal guess = new BigDecimal (initialGuess.toString ());

     // Iterate

     iterations = 0;
     boolean more = true;
     while (more) {
       lastGuess = guess;
       guess = n.divide(guess, scale, BigDecimal.ROUND_HALF_UP);
       guess = guess.add(lastGuess);
       guess = guess.divide (TWO, scale, BigDecimal.ROUND_HALF_UP);
      // System.out.println("Next guess " + guess.toString ());
       error = n.subtract (guess.multiply (guess));
       if (++iterations >= maxIterations) {
         more = false;
       }
       else if (lastGuess.equals (guess)) {
         more = error.abs ().compareTo (ONE) >= 0;
       }
     }
     return guess;

   }
}

