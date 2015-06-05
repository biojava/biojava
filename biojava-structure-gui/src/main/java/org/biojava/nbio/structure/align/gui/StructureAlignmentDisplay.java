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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;
import org.biojava.nbio.structure.align.superimpose.MultipleSuperimposer;
import org.biojava.nbio.structure.align.superimpose.ReferenceSuperimposer;
import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class StructureAlignmentDisplay {
	private static final Logger logger = LoggerFactory.getLogger(StructureAlignmentDisplay.class);
   
   /** Display an AFPChain alignment
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
            
      List<Group> hetatms  = StructureTools.getUnalignedGroups(ca1);
      List<Group> hetatms2 = StructureTools.getUnalignedGroups(ca2);
         
      return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2, hetatms, hetatms2);

   }

   /** 
    * Display a MultipleAlignment. New structures are downloaded if they were not cached in the alignment
    * and they are entirely rotated here with the Pose information.
    * 
    * @param multAln
    * @return MultipleAlignmentJmol instance
    * @throws StructureException
    * @throws StructureAlignmentException 
    * @throws IOException 
    */
   public static MultipleAlignmentJmol display(MultipleAlignment multAln) throws StructureException, StructureAlignmentException, IOException {

	   int size = multAln.size();

	   List<Atom[]> atomArrays = multAln.getEnsemble().getAtomArrays();
	   for (int i=0; i<size; i++){
		   if (atomArrays.get(i).length < 1) 
			   throw new StructureException("Length of atoms arrays is too short! " + atomArrays.get(i).length);
	   }

	   List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();

	   List<Matrix4d> transformations = multAln.getTransformations();
	   if( transformations == null ) {
		   //TODO temporary hack for missing transformations
		   logger.error("BlockSet transformations are unimplemented. Superimposing to first structure.");
		   // clone input, since we're about to re-superimpose it
		   multAln = multAln.clone();
		   MultipleSuperimposer imposer = new ReferenceSuperimposer();
		   imposer.superimpose(multAln);
		   transformations = multAln.getTransformations();
		   assert(transformations != null);
	   }

	   //Rotate the atom coordinates of all the structures
	   for (int i=0; i<size; i++){
		   //TODO handle BlockSet-level transformations for flexible alignments.
		   // In general, make sure this method has the same behavior as the other display. -SB 2015-06

		   // Assume all atoms are from the same structure
		   Structure displayS = atomArrays.get(i)[0].getGroup().getChain().getParent().clone();
		   Atom[] rotCA = StructureTools.getRepresentativeAtomArray(displayS);
		   //Rotate the structure to ensure a full rotation in the display
		   Calc.transform(rotCA[0].getGroup().getChain().getParent(), multAln.getTransformations().get(i));
		   rotatedAtoms.add(rotCA);
	   }


	   MultipleAlignmentJmol jmol = new MultipleAlignmentJmol(multAln, rotatedAtoms);
	   jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
	   return jmol;
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
      
      // List of groups to be rotated according to the alignment
      Group[] twistedGroups = new Group[ ca2.length];
      
      //int blockNum = afpChain.getBlockNum();
            
      int i = -1;
     
      // List of groups from the structure not included in ca2 (e.g. ligands)
      // Will be rotated according to first block
      List<Group> hetatms2 = StructureTools.getUnalignedGroups(ca2);

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

         // Superimpose ligands relative to the first block
         if( hetatms2.size() > 0 ) {
          
            if ( afpChain.getBlockRotationMatrix().length > 0 ) {

               Matrix m1      = afpChain.getBlockRotationMatrix()[0];
               //m1.print(3,3);
               Atom   vector1 = afpChain.getBlockShiftVector()[0];
               //System.out.println("shift vector:" + vector1);

               for ( Group g : hetatms2){                       
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
   public static void shiftCA2(AFPChain afpChain, Atom[] ca2,  Matrix m, Atom shift, Group[] twistedGroups) {
	   
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
