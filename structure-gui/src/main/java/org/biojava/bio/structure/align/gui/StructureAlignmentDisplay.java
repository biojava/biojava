package org.biojava.bio.structure.align.gui;

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.AFPTwister;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

public class StructureAlignmentDisplay {

   public static StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1,
         Atom[] ca2, List<Group> hetatms, List<Group> nucs1,
         List<Group> hetatms2, List<Group> nucs2) throws StructureException {

      Group[] twistedGroups = new Group[ ca2.length];


      if ( afpChain.getBlockRotationMatrix().length == 0 ) {
         throw new StructureException("No rotation matrix found to rotate 2nd structure!");
      }

      int blockNum = afpChain.getBlockNum();
      int i = -1;
     
      if  ( blockNum == 1 ) {
         
         Matrix m   =  afpChain.getBlockRotationMatrix()[ 0];
         Atom shift =  afpChain.getBlockShiftVector()   [ 0 ];

         shiftAll(afpChain, ca2,hetatms2, nucs2, m,shift, twistedGroups);
         
      } else {
         
         for (Atom a: ca2){
            i++;
            twistedGroups[i]=a.getParent();
            
         }


         twistedGroups = AFPTwister.twistOptimized(afpChain, ca1, ca2);


      }


      return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2,hetatms, nucs1, hetatms2, nucs2);

   }

   private static void shiftBlock(AFPChain afpChain, Atom[] ca2, int blockNr, Matrix m, Atom shift, Group[] twistedGroups)
   {

      int[][][] optAln = afpChain.getOptAln();

      int[] optLen = afpChain.getOptLen();



      for(int j = 0; j < optLen[blockNr]; j ++) {
         //p1 = optAln[blockNr][0][j];
         int p2 = optAln[blockNr][1][j];
         Atom a = ca2[p2];
         Calc.rotate(a,m);
         Calc.shift(a,shift);        
      }


   }

   private static void shiftAll(AFPChain afpChain, Atom[] ca2, List<Group> hetatms2, List<Group> nucs2, Matrix m, Atom shift, Group[] twistedGroups)
   {
      int i = -1;
      for (Atom a: ca2){
         i++;
         Group g = a.getParent();
         twistedGroups[i]=g;   
         Calc.rotate(g,m);
         Calc.shift(g, shift);
      }

      //shiftBlock(afpChain, ca2, 0, m, shift, twistedGroups);
//
      // this actually is done by the DisplayAFP ...
//      for (Group g : hetatms2){
//         Calc.rotate(g,m);
//         Calc.shift(g, shift);
//      }
//
//      for (Group g : nucs2){
//         Calc.rotate(g,m);
//         Calc.shift(g, shift);
//      }

   }


}
