package org.biojava.bio.structure.align.model;

import java.io.Serializable;

import org.biojava.bio.structure.jama.Matrix;

/** A class to represent a FATCAT AFP
 *
 * @author Andreas Prlic
 *
 */

public class AFP implements Serializable {

   /**
    * 
    */
   private static final long serialVersionUID = 3901209995477111829L;
   int p1;
   int p2;
   int fragLen;
   double rmsd;
   Matrix m;
   double[] t;
   double score;

   long id;

   public String toString(){

      // we use the metric of

      // Manfred J. Sippl
      // On Distance and Similarity in Fold Space
      // Bioinformatics, 24, pp. 872-873  (2008)


      StringBuffer buf = new StringBuffer();
      buf.append("AFP: p1:");
      buf.append(p1);
      buf.append(" p2: ");
      buf.append(p2);
      buf.append(" len " );
      buf.append(fragLen);
      buf.append(" rmsd ");
      buf.append(rmsd);
      buf.append(" score ");
      buf.append(score);
      return buf.toString();
   }

 
   public long getId()
   {
      return id;
   }

   public void setId(long id)
   {
      this.id = id;
   }

   public int getP1() {
      return p1;
   }
   public void setP1(int p1) {
      this.p1 = p1;
   }
   public int getP2() {
      return p2;
   }
   public void setP2(int p2) {
      this.p2 = p2;
   }
   public int getFragLen() {
      return fragLen;
   }
   public void setFragLen(int fragLen) {
      this.fragLen = fragLen;
   }
   public double getRmsd() {
      return rmsd;
   }
   public void setRmsd(double rmsd) {
      this.rmsd = rmsd;
   }
   public Matrix getM() {
      return m;
   }
   public void setM(Matrix m) {
      this.m = m;
   }
   public double[] getT() {
      return t;
   }
   public void setT(double[] t) {
      this.t = t;
   }
   public double getScore() {
      return score;
   }
   public void setScore(double score) {
      this.score = score;
   }




}
