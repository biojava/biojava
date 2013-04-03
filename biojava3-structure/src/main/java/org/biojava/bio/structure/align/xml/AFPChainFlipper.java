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
 * Created on Sep 9, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.xml;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

public class AFPChainFlipper {


	/** Flip the position of name1 and name2 (as well as all underlying data) in an AFPChain.
	 * This is a utility function for AFPChainXMLParser.
	 * You will have to call AFPCHainXMLParser.rebuildAFPChain in order to get twisted groups...
	 * 
	 * @param o ... the original AFPCHain that should be flipped
	 * @return a cloned AFPCHain which the positions of name1 and name2 flipped.
	 */
	public static AFPChain flipChain(AFPChain o) throws StructureException{

		AFPChain n = new AFPChain();
		n.setAlgorithmName(o.getAlgorithmName());
		n.setVersion(o.getVersion());

		n.setName2(o.getName1());
		n.setName1(o.getName2());

		n.setCa1Length(o.getCa2Length());
		n.setCa2Length(o.getCa1Length());

		int[] optLen = o.getOptLen();
		n.setOptLen(optLen);

		int blockNum = o.getBlockNum();
		n.setBlockNum(blockNum);
		n.setBlockSize(o.getBlockSize());
		n.setBlockScore(o.getBlockScore());
		n.setBlockRmsd(o.getBlockRmsd());
		n.setBlockGap(o.getBlockGap());

		int minLength = Math.min(n.getCa1Length(),n.getCa2Length());
		int[][][] optAlnN 			= new int[blockNum][2][minLength];
		int[][][] optAlnO           = o.getOptAln();


		String[][][] pdbAlnN         = new String[blockNum][2][minLength];
		String[][][] pdbAlnO         = o.getPdbAln();

		if ( ( optAlnO == null) && ( pdbAlnO == null) ){
			System.err.println("Can't get either optAln or pdbAln data from original AFPChain. Not enough information to recreate alignment!");
		}



		for (int blockNr = 0 ; blockNr < blockNum ; blockNr++) {
			for ( int eqrNr = 0 ; eqrNr < optLen[blockNr] ; eqrNr++ ) {

				if ( optAlnO != null ){
					optAlnN[blockNr][0][eqrNr] = optAlnO[blockNr][1][eqrNr];
					optAlnN[blockNr][1][eqrNr] = optAlnO[blockNr][0][eqrNr];
				}
				if ( pdbAlnO != null) {
					pdbAlnN[blockNr][0][eqrNr] = pdbAlnO[blockNr][1][eqrNr];
					pdbAlnN[blockNr][1][eqrNr] = pdbAlnO[blockNr][0][eqrNr];
				}
			}
		}

		n.setOptAln(optAlnN);

		if ( pdbAlnO != null) {
			n.setPdbAln(pdbAlnN);
		}



		n.setAlnLength(o.getAlnLength());
		n.setAlignScore(o.getAlignScore());
		n.setAlignScoreUpdate(o.getAlignScoreUpdate());
		n.setAfpSet(o.getAfpSet());
		n.setChainRmsd(o.getChainRmsd());
		n.setFocusRes1(o.getFocusRes2());
		n.setFocusRes2(o.getFocusRes1());
		n.setFocusResn(o.getFocusResn());
		n.setGapLen(o.getGapLen());
		n.setIdentity(o.getIdentity());
		n.setNormAlignScore(o.getNormAlignScore());
		n.setOptLength(o.getOptLength());
		n.setProbability(o.getProbability());
		n.setSimilarity(o.getSimilarity());
		n.setTotalLenIni(o.getTotalLenIni());		
		n.setTotalRmsdIni(o.getTotalRmsdIni());
		n.setTotalRmsdOpt(o.getTotalRmsdOpt());
		n.setTMScore(o.getTMScore());
		
		
		// change direction of the Matrix and shift!
		// 
		Matrix[] maxO  = o.getBlockRotationMatrix();
		Matrix[] maxN = new Matrix[maxO.length];

		int i = -1;

		Atom[] shiftO = o.getBlockShiftVector();
		Atom[] shiftN = new Atom[shiftO.length];

		for (Matrix m : maxO){
			i++;
			if ( m == null) {
				// alignment too short probably
				continue;
			}
			try {
				Matrix mnew = m ;
				Atom a = shiftO[i];
				
				maxN[i] = mnew.transpose();

				shiftN[i] =  Calc.invert(a);
				
				 Calc.rotate(shiftN[i],maxN[i]);

			} catch (StructureException e){
				e.printStackTrace();
			}
		}

		n.setBlockRotationMatrix(maxN);
		n.setBlockShiftVector(shiftN);
		return n;

	}
}