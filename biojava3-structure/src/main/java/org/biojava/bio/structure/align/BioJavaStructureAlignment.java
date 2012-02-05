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
 * Created on Oct 2, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

import org.biojava.bio.structure.align.helper.IndexPair;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.jama.Matrix;


/** Wrapper for the BioJava Structure Alignment Implementation
 * 
 * @author Andreas Prlic
 *
 */
public class BioJavaStructureAlignment

implements StructureAlignment  {

	public static final String algorithmName = "BioJava_structure";
	private static final float versionNr = 0.1f;
	StrucAligParameters params;

	public BioJavaStructureAlignment(){
		params = new StrucAligParameters();
	}

//	public StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1,
//			Atom[] ca2, List<Group> hetatms, List<Group> nucs1,
//			List<Group> hetatms2, List<Group> nucs2) throws StructureException {
//		
//	
//
//		Group[] twistedGroups = new Group[ ca2.length];
//		
//		int i=-1;
//		Matrix m =  afpChain.getBlockRotationMatrix()[0];
//		Atom shift =  afpChain.getBlockShiftVector()[0];
//		
//		for (Atom a: ca2){
//			i++;
//			twistedGroups[i]=a.getParent();
//			Calc.rotate(twistedGroups[i],m);
//			Calc.shift(twistedGroups[i],shift);
//		}
//		
//		return DisplayAFP.display(afpChain, twistedGroups, ca1, ca2,hetatms, nucs2, hetatms2, nucs2);
//	}

	public String getAlgorithmName() {
		return algorithmName;
	}

	public ConfigStrucAligParams getParameters() {
		return null;//TODO shall we update it?
	}
	
	public void setParameters(ConfigStrucAligParams o){
		//TODO what is the relation between StrucAligParameters and ConfigStrucAligParams?
		if ( ! (o instanceof StrucAligParameters)){
			throw new IllegalArgumentException("Provided parameters are not of type StrucAligParameters!");
		}
		params = (StrucAligParameters) o;
	}

	public String getVersion() {
		return versionNr+"";
	}

	public String printHelp() {
		return "not implemented yet. Algorithm still under development.";
	}



	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
		StrucAligParameters params = StrucAligParameters.getDefaultParameters();
		return align(ca1,ca2,params);
		
	}


	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
	throws StructureException {
		if ( ! (params instanceof StrucAligParameters)){
			throw new IllegalArgumentException("BioJava structure alignment requires a StrucAligParameters class for the arguments.");
		}
		this.params = (StrucAligParameters) params;

		AFPChain afpChain = new AFPChain();
		StructurePairAligner aligner = new StructurePairAligner();
		aligner.align(ca1,ca2,this.params);

		// copy the results over into the AFPChain...
		AlternativeAlignment[] aligs = aligner.getAlignments();
		if ( aligs.length > 0){

			AlternativeAlignment altAlig = aligs[0];
			// copy the results over!
			copyResults(afpChain,altAlig,ca1,ca2);


		}

		return afpChain;

	}



	private void copyResults(AFPChain afpChain, AlternativeAlignment altAlig, Atom[] ca1, Atom[] ca2) {
		afpChain.setAlgorithmName(getAlgorithmName());
		afpChain.setVersion(getVersion());
		afpChain.setAlignScore(altAlig.getScore());
		afpChain.setOptLength(altAlig.getEqr());
		afpChain.setBlockRotationMatrix(new Matrix[]{altAlig.getRotationMatrix()});
		afpChain.setBlockShiftVector(new Atom[]{altAlig.getShift() });
		afpChain.setBlockNum(1);
		
		double rmsd = altAlig.getRmsd();
		afpChain.setBlockRmsd(new double[]{rmsd});
		
		
		int nAtom = altAlig.getEqr();
		int lcmp = altAlig.getPath().length;
		int[] optLen = new int[]{nAtom};
		afpChain.setOptLen(optLen);
		afpChain.setOptLength(nAtom);		
		afpChain.setAlnLength(lcmp);
				
		int[][][] optAln = new int[1][2][lcmp];
		afpChain.setOptAln(optAln);
		
		afpChain.setOptRmsd(new double[]{rmsd});
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setChainRmsd(rmsd);
		//afpChain.setProbability(-1);
		int nse1 = ca1.length;
		int nse2 = ca2.length;
		
		char[] alnseq1 = new char[nse1+nse2+1];
		char[] alnseq2 = new char[nse1+nse2+1] ;
		char[] alnsymb = new char[nse1+nse2+1];
     
		IndexPair[] path = altAlig.getPath(); 
        
		int pos = 0;
		for(int ia=0; ia<lcmp; ia++) {

			IndexPair align_se = path[ia];
			// no gap
			if(align_se.getRow() !=-1 && align_se.getCol()!=-1) {

				optAln[0][0][pos] = align_se.getRow();
				optAln[0][1][pos] = align_se.getCol();
				
				char l1 = getOneLetter(ca1[align_se.getRow()].getGroup());
				char l2 = getOneLetter(ca2[align_se.getCol()].getGroup());
				
				alnseq1[ia] = Character.toUpperCase(l1);
				alnseq2[ia] = Character.toUpperCase(l2);
				alnsymb[ia] = '1';
				pos++;
								
			} else {
				// there is a gap on this position
				alnsymb[ia] = ' ';
				if (align_se.getRow() == -1 ) {
					alnseq1[ia] = '-';
				} else {
					char l1 = getOneLetter(ca1[align_se.getRow()].getGroup());
					alnseq1[ia] = Character.toLowerCase(l1);
				}
				if ( align_se.getCol() == -1 ) {
					alnseq2[ia] = '-';
				} else {
					char l2 = getOneLetter(ca2[align_se.getCol()].getGroup());
					alnseq2[ia] = Character.toLowerCase(l2);
				}
				
			}
		}
		afpChain.setAlnseq1(alnseq1);
        afpChain.setAlnseq2(alnseq2);
        afpChain.setAlnsymb(alnsymb);

	}
	
	private static char getOneLetter(Group g){

        try {
           Character c = StructureTools.get1LetterCode(g.getPDBName());
           return c;
        } catch (Exception e){
           return 'X';
        }
     }

}
