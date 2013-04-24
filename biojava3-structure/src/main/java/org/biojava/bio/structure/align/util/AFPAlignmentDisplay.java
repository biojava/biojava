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
 * Created on Jun 17, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.GuiWrapper;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.jama.Matrix;


public class AFPAlignmentDisplay
{

	private static final int[][] aaMatrix = new int[][] {{6,0,-2,-3,-2,0,-1,0,-2,-2,-2,-2,-2,-3,-4,-4,-3,-3,-3,-2,-4},
		{0,4,-1,0,0,1,-2,-2,-1,-1,-1,-2,-1,0,-1,-1,-1,-2,-2,-3,-4},
		{-2,-1,7,-3,-1,-1,-1,-2,-1,-1,-1,-2,-2,-2,-3,-3,-2,-4,-3,-4,-4},
		{-3,0,-3,9,-1,-1,-3,-3,-4,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2,-4},
		{-2,0,-1,-1,5,1,-1,0,-1,-1,-1,-2,-1,0,-1,-1,-1,-2,-2,-2,-4},
		{0,1,-1,-1,1,4,0,1,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-3,-4},
		{-1,-2,-1,-3,-1,0,6,1,2,0,-1,-1,-2,-3,-3,-4,-3,-3,-3,-4,-4},
		{0,-2,-2,-3,0,1,1,6,0,0,0,1,0,-3,-3,-3,-2,-3,-2,-4,-4},
		{-2,-1,-1,-4,-1,0,2,0,5,2,1,0,0,-2,-3,-3,-2,-3,-2,-3,-4},
		{-2,-1,-1,-3,-1,0,0,0,2,5,1,0,1,-2,-3,-2,0,-3,-1,-2,-4},
		{-2,-1,-1,-3,-1,0,-1,0,1,1,5,-1,2,-2,-3,-2,-1,-3,-2,-3,-4},
		{-2,-2,-2,-3,-2,-1,-1,1,0,0,-1,8,0,-3,-3,-3,-2,-1,2,-2,-4},
		{-2,-1,-2,-3,-1,-1,-2,0,0,1,2,0,5,-3,-3,-2,-1,-3,-2,-3,-4},
		{-3,0,-2,-1,0,-2,-3,-3,-2,-2,-2,-3,-3,4,3,1,1,-1,-1,-3,-4},
		{-4,-1,-3,-1,-1,-2,-3,-3,-3,-3,-3,-3,-3,3,4,2,1,0,-1,-3,-4},
		{-4,-1,-3,-1,-1,-2,-4,-3,-3,-2,-2,-3,-2,1,2,4,2,0,-1,-2,-4},
		{-3,-1,-2,-1,-1,-1,-3,-2,-2,0,-1,-2,-1,1,1,2,5,0,-1,-1,-4},
		{-3,-2,-4,-2,-2,-2,-3,-3,-3,-3,-3,-1,-3,-1,0,0,0,6,3,1,-4},
		{-3,-2,-3,-2,-2,-2,-3,-2,-2,-1,-2,2,-2,-1,-1,-1,-1,3,7,2,-4},
		{-2,-3,-4,-2,-2,-3,-4,-4,-3,-2,-3,-2,-3,-3,-3,-2,-1,1,2,11,-4},
		{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1}};

	private static Character[] aa1 = new Character[]{  'G',   'A',   'P',   'C',   'T',   'S','D',   'N',   'E',   'Q',   'K',   'H',   'R', 'V',   'I',   'L',   'M',   'F',   'Y',   'W', '-'};

	private static final List<Character> aa1List = Arrays.asList(aa1);

	public static Matrix getRotMax(AFPChain afpChain,Atom[] ca1,Atom[] ca2) throws StructureException{

		Atom[] a1 = getAlignedAtoms1(afpChain,ca1);
		Atom[] a2 = getAlignedAtoms2(afpChain,ca2);

		SVDSuperimposer svd = new SVDSuperimposer(a1,a2);

		return svd.getRotation();

	}

	public static Atom getTranslation(AFPChain afpChain,Atom[] ca1,Atom[] ca2) throws StructureException{


		Atom[] a1 = getAlignedAtoms1(afpChain,ca1);
		Atom[] a2 = getAlignedAtoms2(afpChain,ca2);

		SVDSuperimposer svd = new SVDSuperimposer(a1,a2);

		return svd.getTranslation();

	}

	public static Atom[] getAlignedAtoms1(AFPChain afpChain,Atom[] ca1){
		List<Atom> atoms = new ArrayList<Atom>();

		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();


		for(int bk = 0; bk < blockNum; bk ++)       {


			for ( int i=0;i< optLen[bk];i++){
				int pos = optAln[bk][0][i];
				atoms.add(ca1[pos]);
			}

		}
		return (Atom[])atoms.toArray(new Atom[atoms.size()]);
	}
	public static Atom[] getAlignedAtoms2(AFPChain afpChain,Atom[] ca2){

		List<Atom> atoms = new ArrayList<Atom>();

		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();


		for(int bk = 0; bk < blockNum; bk ++)       {


			for ( int i=0;i< optLen[bk];i++){
				int pos = optAln[bk][1][i];
				atoms.add(ca2[pos]);
			}

		}
		return (Atom[])atoms.toArray(new Atom[atoms.size()]);
	}



	/**
     * Extract the alignment output
     * <p>eg<pre>
     * STWNTWACTWHLKQP--WSTILILA
     * 111111111111     22222222
     * SQNNTYACSWKLKSWNNNSTILILG
     * </pre>
     * Those position pairs labeled by 1 and 2 are equivalent positions, belongs to
     * two blocks 1 and 2. The residues between labeled residues are non-equivalent,
     * with '-' filling in their resulting gaps.
     * <p>
     * The three lines can be accessed using 
     * {@link AFPChain#getAlnseq1()}, {@link AFPChain#getAlnsymb()},
     * and {@link AFPChain#getAlnseq2()}.
     * 
	 */
	public static void getAlign(AFPChain afpChain,Atom[] ca1,Atom[] ca2)
	{
		boolean showSeq = false;
		getAlign(afpChain, ca1, ca2, showSeq);
	}

	/**
	 * Sets the following properties:
	 * <ul>
	 * <li>The alignment strings {@link AFPChain#setAlnseq1(char[]) alnseq1},
	 *  {@link AFPChain#setAlnseq2(char[]) alnseq2},
	 *  and {@link AFPChain#setAlnsymb(char[]) alnsymb}</li>
	 * <li>{@link AFPChain#setAlnbeg1(int) alnbeg1} and 2</li>
	 * <li>{@link AFPChain#setAlnLength(int) alnLength} and 
	 *  {@link AFPChain#setGapLen(int) gapLen}</li>
	 * </ul>
	 * <p>
	 * Expects the following properties to be previously computed:
	 * <ul>
	 * <li>{@link AFPChain#getOptAln()} and lengths
	 * </ul>
	 * 
	 * <section>Known Bugs</section>
	 * Expects the alignment to have linear topology. May give odd results
	 * for circular permutations and other complicated topologies.
	 * 
	 * @param afpChain Alignment between ca1 and ca2
	 * @param ca1 CA atoms of the first protein
	 * @param ca2 CA atoms of the second protein
	 * @param showSeq Use symbols reflecting sequence similarity: '|' for identical,
	 *  ':' for similar, '.' for dissimilar. Otherwise, use the block number
	 *  to show aligned residues.
	 */
	public static void getAlign(AFPChain afpChain,Atom[] ca1,Atom[] ca2, boolean showSeq) {

		char[] alnsymb = afpChain.getAlnsymb();
		char[] alnseq1 = afpChain.getAlnseq1();
		char[] alnseq2 = afpChain.getAlnseq2();

		int     i, j, k, p1, p2, p1b, p2b, lmax;

		int pro1Len = ca1.length;
		int pro2Len = ca2.length;

		p1b = p2b = 0;



		if(alnsymb == null)     {
			alnseq1 = new char[pro1Len+pro2Len +1];
			alnseq2 = new char[pro1Len+pro2Len +1] ;
			alnsymb = new char[pro1Len+pro2Len+1];

			afpChain.setAlnseq1(alnseq1);
			afpChain.setAlnseq2(alnseq2);
			afpChain.setAlnsymb(alnsymb);
		}



		//		if ( afpChain.getAlgorithmName().equals(FatCatRigid.algorithmName)){
		//			showSeq = true;
		//		}

		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		int alnbeg1 = afpChain.getAlnbeg1(); // immediately overwritten
		int alnbeg2 = afpChain.getAlnbeg2(); // immediately overwritten
		int alnLength = afpChain.getAlnLength(); // immediately overwritten
		int optLength = afpChain.getOptLength();

		if ( optLen == null) {
			optLen = new int[blockNum];
			for (int oi =0 ; oi < blockNum ; oi++)
				optLen[oi] = 0;
		}
		int     len = 0;
		for(i = 0; i < blockNum; i ++)  {        	
			for(j = 0; j < optLen[i]; j ++) {
				p1 = optAln[i][0][j];
				p2 = optAln[i][1][j];

				// weird, could not find a residue in the Atom array. Did something change in the underlying data?
				if (( p1 == -1 ) || ( p2 == -1)) {
					System.err.println("Could not get atom on position " + j );
					continue;
				}
				if(len > 0)     {
					lmax = (p1 - p1b - 1)>(p2 - p2b - 1)?(p1 - p1b - 1):(p2 - p2b - 1);
					for(k = 0; k < lmax; k ++)      {
						if(k >= (p1 - p1b - 1)) alnseq1[len] = '-';
						else {
							char oneletter = getOneLetter(ca1[p1b+1+k].getGroup());
							alnseq1[len] = oneletter;
						}
						if(k >= (p2 - p2b - 1)) alnseq2[len] = '-';
						else  {
							char oneletter = getOneLetter(ca2[p2b+1+k].getGroup());
							alnseq2[len] = oneletter;
						}
						alnsymb[len ++] = ' ';
					}
				}
				else    {
					alnbeg1 = p1; //the first position of sequence in alignment
					alnbeg2 = p2;
				}
				
				if ( p1 < ca1.length && p2<ca2.length) {
				
					alnseq1[len] = getOneLetter(ca1[p1].getGroup());              
					alnseq2[len] = getOneLetter(ca2[p2].getGroup());
				} else {
					//TODO handle permutations
					alnseq1[len]='?';
					alnseq2[len]='?';
				}
				if ( showSeq) {
					if ( alnseq1[len] == alnseq2[len]){
						alnsymb[len ++] = '|';
					} else {
						double score = aaScore(alnseq1[len], alnseq2[len]);

						if ( score > 1)
							alnsymb[len ++] = ':';
						else 
							alnsymb[len ++] = '.';
					}
				} else {
					String tmpS = String.format( "%d", i + 1);
					alnsymb[len ++] = tmpS.charAt(0);
				}
				p1b = p1;
				p2b = p2;
			}
		}
		alnLength = len;


		//		System.out.println(alnseq1);
		//		System.out.println(alnsymb);
		//		System.out.println(alnseq2);

		afpChain.setOptAln(optAln);
		afpChain.setOptLen(optLen);
		afpChain.setAlnbeg1(alnbeg1);
		afpChain.setAlnbeg2(alnbeg2);
		afpChain.setAlnLength(alnLength);
		afpChain.setGapLen(alnLength-optLength);
	}

	private static char getOneLetter(Group g){

		try {
			Character c = StructureTools.get1LetterCode(g.getPDBName());
			return c;
		} catch (Exception e){
			return 'X';
		}
	}

	public static int aaScore(char a, char b){
		if (a == 'x')
			a= '-';
		if (b == 'x')
			b='-';
		if ( a == 'X')
			a = '-';
		if ( b == 'X')
			b = '-';

		int pos1 = aa1List.indexOf(a);
		int pos2 = aa1List.indexOf(b);


		// SEC an PYL amino acids are not supported as of yet...

		if ( pos1 < 0) {
			System.err.println("unknown char " + a);
			return 0;
		}
		if ( pos2 < 0) {
			System.err.println("unknown char " + b);
			return 0;
		}

		return aaMatrix[pos1][pos2];
	}

	public static Map<String,Double> calcIdSimilarity(char[] seq1, char[] seq2, int alnLength){
		double identity = 0.0;
		double similarity = 0.0;
		
		if ( seq1 == null || seq2 == null){
			System.err.println("AFPAlignmentDisplay can't calc ID if alignment strings are null! ");
			Map<String, Double> m = new HashMap<String, Double>();
			m.put("similarity", similarity);
			m.put("identity", identity);
			return m;
		}
		
		int     i;

		for(i = 0; i < alnLength; i ++)        {
		
				if(seq1[i] == seq2[i])  {

					identity += 1.0;
				}
				else if(seq1[i] == '-' || seq1[i] == '*' || seq1[i] == '.' ||
						seq2[i] == '-' || seq2[i] == '*' || seq2[i] == '.' )
					continue;
				else if(aaScore(seq1[i], seq2[i]) > 0)  similarity += 1.0;		
		}


		if ( alnLength > 0){
			similarity = (identity + similarity) / (double)alnLength;
			identity = identity/(double)alnLength;
		} else {
			
		}
		Map<String, Double> m = new HashMap<String, Double>();
		m.put("similarity", similarity);
		m.put("identity", identity);
		
		return m;
	}


	public static Structure createArtificalStructure(AFPChain afpChain, Atom[] ca1,
			Atom[] ca2) throws Exception{


		if ( afpChain.getNrEQR() < 1){
			return GuiWrapper.getAlignedStructure(ca1, ca2);
		}

		Group[] twistedGroups = GuiWrapper.prepareGroupsForDisplay(afpChain,ca1, ca2);

		List<Atom> twistedAs = new ArrayList<Atom>();

		for ( Group g: twistedGroups){
			if ( g == null )
				continue;
			if ( g.size() < 1)
				continue;
			Atom a = g.getAtom(0);
			twistedAs.add(a);
		}
		Atom[] twistedAtoms = (Atom[])twistedAs.toArray(new Atom[twistedAs.size()]);

		List<Group> hetatms  = new ArrayList<Group>();
		List<Group> nucs1    = new ArrayList<Group>();
		Group g1 = ca1[0].getGroup();
		Chain c1 = null;
		if ( g1 != null) {
			c1 = g1.getChain();
			if ( c1 != null){
				hetatms = c1.getAtomGroups("hetatm");;
				nucs1  = c1.getAtomGroups("nucleotide");
			}
		}
		List<Group> hetatms2 = new ArrayList<Group>();
		List<Group> nucs2    = new ArrayList<Group>();
		Group g2 = ca2[0].getGroup();
		Chain c2 = null;
		if ( g2 != null){
			c2 = g2.getChain();
			if ( c2 != null){
				hetatms2 = c2.getAtomGroups("hetatm");
				nucs2 = c2.getAtomGroups("nucleotide");
			}
		}
		Atom[] arr1 = GuiWrapper.getAtomArray(ca1, hetatms, nucs1);
		Atom[] arr2 = GuiWrapper.getAtomArray(twistedAtoms, hetatms2, nucs2);

		Structure artificial = GuiWrapper.getAlignedStructure(arr1,arr2);
		return artificial;
	}
	
	/** get the block number for an aligned position
	 * 
	 * @param afpChain
	 * @param aligPos
	 * @return
	 */
	public static int getBlockNrForAlignPos(AFPChain afpChain, int aligPos){
		// moved here from DisplayAFP;
		
		//int ungappedPos = -1;
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();

		int len = 0;
		int p1b=0;
		int p2b=0;

		for(int i = 0; i < blockNum; i ++)  {   

			for(int j = 0; j < optLen[i]; j ++) {

				int p1 = optAln[i][0][j];
				int p2 = optAln[i][1][j];

				//                 System.out.println(p1 + " " + p2 + " " +  footer2.toString());

				if ( len == 0){
					//
				} else {
					// check for gapped region
					int lmax = (p1 - p1b - 1)>(p2 - p2b - 1)?(p1 - p1b - 1):(p2 - p2b - 1);
					for(int k = 0; k < lmax; k ++)      {
						len++;
					}
				}

				len++;
				//ungappedPos++;
				p1b = p1;
				p2b = p2;
				if ( len >= aligPos) {

					return i;
				}
			}
		}

		return blockNum;

	}


}
