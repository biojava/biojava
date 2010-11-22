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
 * Created on May 24, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.bio.structure.align.pairwise;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/** A class to track the alignment results in a flat file
 *
 * @author andreas
 *
 */
public class AlignmentResult implements Serializable
{

   /**
    *
    */
   private static final long serialVersionUID = -4132105905712445473L;
   AlternativeAlignment[] alignments;
   String pdb1;
   String pdb2;

   String chain1;
   String chain2;

   int length1;
   int length2;

   long calculationTime;
   long ioTime;

   public String toString(){
      StringBuffer buf = new StringBuffer();
      buf.append(pdb1);
      buf.append("_");
      buf.append(chain1);
      buf.append(" vs. ");
      buf.append(pdb2);
      buf.append("_");
      buf.append(chain2);
      buf.append(" : ");
      buf.append(" l1: ");
      buf.append(length1);
      buf.append(" l2: ");
      buf.append(length2);
      buf.append(" ");
      if ( alignments != null)
         if ( alignments.length > 0) {
            AlternativeAlignment a = alignments[0];
            buf.append(a.toString());
            int eqr = a.getEqr();
            buf.append(" %res1: ");
            buf.append(Math.round((eqr/(float)length2)*100));
            buf.append(" %res2: ");
            buf.append(Math.round((eqr/(float)length1)*100));
            buf.append(" ");

         }
      buf.append(" ioTime: ");
      buf.append(ioTime);
      buf.append(" compTime: ");
      buf.append(calculationTime);
      return buf.toString();

   }

   public AlternativeAlignment[] getAlignments()
   {
      return alignments;
   }

   /** we only keep the first alternative...
    *
    * @param alignments
    */
   public void setAlignments(AlternativeAlignment[] alignments)
   {
      if ( alignments.length > 0){
         this.alignments = new AlternativeAlignment[1];
         this.alignments[0]=alignments[0];
      }

   }
   public String getPdb1()
   {
      return pdb1;
   }
   public void setPdb1(String pdb1)
   {
      this.pdb1 = pdb1;
   }
   public String getPdb2()
   {
      return pdb2;
   }
   public void setPdb2(String pdb2)
   {
      this.pdb2 = pdb2;
   }
   public String getChain1()
   {
      return chain1;
   }
   public void setChain1(String chain1)
   {
      this.chain1 = chain1;
   }
   public String getChain2()
   {
      return chain2;
   }
   public void setChain2(String chain2)
   {
      this.chain2 = chain2;
   }
   public int getLength1()
   {
      return length1;
   }
   public void setLength1(int length1)
   {
      this.length1 = length1;
   }
   public int getLength2()
   {
      return length2;
   }
   public void setLength2(int length2)
   {
      this.length2 = length2;
   }


   public long getCalculationTime()
   {
      return calculationTime;
   }
   public void setCalculationTime(long calculationTime)
   {
      this.calculationTime = calculationTime;
   }
   public long getIoTime()
   {
      return ioTime;
   }
   public void setIoTime(long ioTime)
   {
      this.ioTime = ioTime;
   }
   public void serialize (File output)
   throws FileNotFoundException, IOException{
      // save alignment result:

      FileOutputStream outStream = new FileOutputStream(output);
      ObjectOutputStream objStream = new ObjectOutputStream(outStream);
      objStream.writeObject(this);
      objStream.close();

   }

   public static AlignmentResult deserialize(File output)
   throws FileNotFoundException, IOException, ClassNotFoundException{
      FileInputStream fin = new FileInputStream(output);
      ObjectInputStream objIn = new ObjectInputStream(fin);
      AlignmentResult result = (AlignmentResult) objIn.readObject();
      objIn.close();
      return result;
   }


}
