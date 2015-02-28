/*
 *                    BioJava development code
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
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.jama.Matrix;

import java.util.ArrayList;
import java.util.List;

public class StructureAlignmentDisplay {

   
   /** Display the alignment
    * 
    * @param afpChain
    * @param ca1
    * @param ca2
    * @return a StructureAlignmentJmol instance
    * @throws StructureException
    */
   public static StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException {
      
      if ( ca1.length < 1 || ca2.length < 1){
         throw new StructureException("length of atoms arrays is too short! " + ca1.length + "," + ca2.length);
      }
      
      Group[] twistedGroups = prepareGroupsForDisplay(afpChain, ca1, ca2);
            
      List<Group> hetatms  = new ArrayList<Group>();
      List<Group> nucs1    = new ArrayList<Group>();
      Group g1 = ca1[0].getGroup();
      Chain c1 = null;
      if ( g1 != null) {
         c1 = g1.getChain();
         if ( c1 != null){
            hetatms = c1.getAtomGroups(GroupType.HETATM);;
            nucs1  = c1.getAtomGroups(GroupType.NUCLEOTIDE);
         }
      }
      List<Group> hetatms2 = new ArrayList<Group>();
      List<Group> nucs2    = new ArrayList<Group>();
      Group g2 = ca2[0].getGroup();
      Chain c2 = null;
      if ( g2 != null){
         c2 = g2.getChain();
         if ( c2 != null){
            hetatms2 = c2.getAtomGroups(GroupType.HETATM);
            nucs2 = c2.getAtomGroups(GroupType.NUCLEOTIDE);
         }
      }
         
      return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2,hetatms, nucs1, hetatms2, nucs2);

   }
   
   /** Rotate the Atoms/Groups so they are aligned for the 3D visualisation
    * 
    * @param afpChain
    * @param ca1
    * @param ca2
    * @return an array of Groups that are transformed for 3D display
    * @throws StructureException
    */
   public static Group[] prepareGroupsForDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
      

      if ( afpChain.getBlockRotationMatrix().length == 0 ) {
         // probably the alignment is too short!
         System.err.println("No rotation matrix found to rotate 2nd structure!");
         afpChain.setBlockRotationMatrix(new Matrix[]{Matrix.identity(3, 3)});
         afpChain.setBlockShiftVector(new Atom[]{new AtomImpl()});
      }
      
      Group[] twistedGroups = new Group[ ca2.length];
      
      //int blockNum = afpChain.getBlockNum();
            
      int i = -1;
     
   
      List<Group> hetatms2 = new ArrayList<Group>();
      List<Group> nucs2    = new ArrayList<Group>();
      
     
      Group g2 = ca2[0].getGroup();
      Chain c2 = null;
      if ( g2 != null){
         c2 = g2.getChain();
         if ( c2 != null){
            hetatms2 = c2.getAtomGroups(GroupType.HETATM);
            nucs2 = c2.getAtomGroups(GroupType.NUCLEOTIDE);
         }
      }
      
      if (  (afpChain.getAlgorithmName().equals(FatCatRigid.algorithmName) ) || (afpChain.getAlgorithmName().equals(FatCatFlexible.algorithmName) ) ){
         
         for (Atom a: ca2){
            i++;
            twistedGroups[i]=a.getGroup();
            
         }

         twistedGroups = AFPTwister.twistOptimized(afpChain, ca1, ca2);
         
      //} else  if  (( blockNum == 1 ) || (afpChain.getAlgorithmName().equals(CeCPMain.algorithmName))) {
      } else {
         
         Matrix m   =  afpChain.getBlockRotationMatrix()[ 0];
         Atom shift =  afpChain.getBlockShiftVector()   [ 0 ];

         shiftCA2(afpChain, ca2, m,shift, twistedGroups);
       
      }
      
      if ( afpChain.getBlockNum() > 0){

         if (( hetatms2.size() > 0) || (nucs2.size() >0)) {
          
            if ( afpChain.getBlockRotationMatrix().length > 0 ) {

               Matrix m1      = afpChain.getBlockRotationMatrix()[0];
               //m1.print(3,3);
               Atom   vector1 = afpChain.getBlockShiftVector()[0];
               //System.out.println("shift vector:" + vector1);

               for ( Group g : hetatms2){                       
                  Calc.rotate(g, m1);
                  Calc.shift(g,vector1);
               }
               for (Group g: nucs2){
                  Calc.rotate(g, m1);
                  Calc.shift(g,vector1);
               }
            }
         }
      }
      
      return twistedGroups;
   }
   

     

   

  /** only shift CA positions.
   * 
  
   */
   public static void shiftCA2(AFPChain afpChain, Atom[] ca2,  Matrix m, Atom shift, Group[] twistedGroups)
   {
      int i = -1;
      for (Atom a: ca2){
         i++;
         Group g = a.getGroup();
        
         
         Calc.rotate(g,m);
         Calc.shift(g, shift);
         
         if (g.hasAltLoc()){
        	 for (Group alt: g.getAltLocs()){
        		 for (Atom alta : alt.getAtoms()){
        			 if ( g.getAtoms().contains(alta))
        				 continue;
        			 Calc.rotate(alta,m);
        			 Calc.shift(alta,shift);
        		 }
        	 }
         }
         
         
         twistedGroups[i]=g;
      }

    
   }


}
