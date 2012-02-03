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
 * Created on Apr 7, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.util;



import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;

public class AFPChainScorer
{
	
	
	public  static double getTMScore(AFPChain align, Atom[] ca1, Atom[] ca2) throws StructureException
	{
		if ( align.getNrEQR() == 0)
			return -1;

		
		// Create new arrays for the subset of atoms in the alignment.
		Atom[] ca1aligned = new Atom[align.getOptLength()];
		Atom[] ca2aligned = new Atom[align.getOptLength()];
		int pos=0;
		int[] blockLens = align.getOptLen();
		int[][][] optAln = align.getOptAln();
		assert(align.getBlockNum() <= optAln.length);
		
		for(int block=0;block< align.getBlockNum();block++) {
		
			if ( ! ( blockLens[block] <= optAln[block][0].length)) {
				System.err.println("AFPChainScorer getTMScore: errors reconstructing alignment block ["+ block +"]. Length is " + blockLens[block] + " but should be <=" + optAln[block][0].length );
			}

			for(int i=0;i<blockLens[block];i++) {
				int pos1 = optAln[block][0][i];
				int pos2 = optAln[block][1][i];
//				if ( pos1 >= ca1.length) {
//					System.err.println("getTMScore: pos1 " + pos1 + " > ca1.length " + ca1.length + " " + align.getName1() + ":" +align.getName2());
//					System.err.println(ca1[0].getGroup().getChain().getChainID());
//					System.err.println(ca1[0].getGroup().getChain().getParent().getPDBCode());
//					
//				}
//				if ( pos2 >= ca2.length) {
//					System.err.println("getTMScore: pos2 " + pos2 + " > ca2.length " + ca2.length + " " + align.getName2() + ":" + align.getName1());
//					System.err.println(ca2[0].getGroup().getChain().getChainID());
//					System.err.println(ca2[0].getGroup().getChain().getParent().getPDBCode());
//				}
				Atom a1 = ca1[pos1];
				Atom a2 = (Atom) ca2[pos2].clone();
								
				ca1aligned[pos] = a1;
				ca2aligned[pos] = a2;
				pos++;
			}
		}
		
		// this can happen when we load an old XML serialization which did not support modern ChemComp representation of modified residues.		
		if ( pos != align.getOptLength()){
			System.err.println("AFPChainScorer getTMScore: Problems reconstructing alignment! nr of loaded atoms is " + pos + " but should be " + align.getOptLength());
			// we need to resize the array, because we allocated too many atoms earlier on.
			ca1aligned = (Atom[]) resizeArray(ca1aligned, pos);
			ca2aligned = (Atom[]) resizeArray(ca2aligned, pos);
		}
		//Superimpose
		SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
		Matrix matrix = svd.getRotation();
		Atom shift = svd.getTranslation();

		for(Atom a : ca2aligned) {
			Calc.rotate(a, matrix);
			Calc.shift(a, shift);
		}

		return SVDSuperimposer.getTMScore(ca1aligned, ca2aligned, ca1.length, ca2.length);

	}

	/**
	 * Reallocates an array with a new size, and copies the contents
	 * of the old array to the new array.
	 * @param oldArray  the old array, to be reallocated.
	 * @param newSize   the new array size.
	 * @return          A new array with the same contents.
	 */
	private static Object resizeArray (Object oldArray, int newSize) {
		int oldSize = java.lang.reflect.Array.getLength(oldArray);
		@SuppressWarnings("rawtypes")
		Class elementType = oldArray.getClass().getComponentType();
		Object newArray = java.lang.reflect.Array.newInstance(
				elementType,newSize);
		int preserveLength = Math.min(oldSize,newSize);
		if (preserveLength > 0)
			System.arraycopy (oldArray,0,newArray,0,preserveLength);
		return newArray; }

}
