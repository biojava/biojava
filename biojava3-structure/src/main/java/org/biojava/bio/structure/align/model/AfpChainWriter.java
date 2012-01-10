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
 * Created on Feb 15, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.model;

import java.io.StringWriter;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;

/** A class to convert the data in an AfpChain object to various String outputs.
 *  
 * @author Andreas Prlic
 *
 */
public class AfpChainWriter
{

	public static final String newline = System.getProperty("line.separator");

	private static int LINELENGTH = 70;

	public static String toFatCat(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
	{

		boolean printLegend = true;
		boolean longHeader  = true;
		boolean showHTML = false;
		boolean showAlignmentBlock = false;
		
		return toFatCatCore(afpChain, ca1, ca2, printLegend, longHeader, showHTML, showAlignmentBlock);   
	}

	public static String toScoresList(AFPChain afpChain){

		// see sippl On distance and similarity in fold space 2008 bioinformatics

		StringWriter writer = new StringWriter();

		if ( afpChain.getAlgorithmName().startsWith("CE")) {
			writer.append("Z-score " );
			writer.append(String.format("%.2f", afpChain.getProbability()));
			writer.append(newline);
		}


		writer.append("Sab (nr. equivalent residues): " );
        writer.append(String.valueOf(afpChain.getNrEQR())).append("");
		writer.append(newline);

		writer.append("Dab (distance between folds a,b): ");
		int dab = afpChain.getCa1Length()+afpChain.getCa2Length() - 2 * afpChain.getNrEQR();
        writer.append(String.valueOf(dab)).append("");
		writer.append(newline);

		writer.append("sab (relative similarity): ");
		double sab = 2 * afpChain.getNrEQR() / (double)( afpChain.getCa1Length() + afpChain.getCa2Length());
        writer.append(String.valueOf(sab)).append("");
		writer.append(newline);

		writer.append("cab (coverage a): ");
		double cab = afpChain.getNrEQR() / (double) afpChain.getCa1Length();
        writer.append(String.valueOf(cab)).append("");
		writer.append(newline);

		writer.append("cba (coverage b): ");
		double cba = afpChain.getNrEQR() / (double) afpChain.getCa2Length();
        writer.append(String.valueOf(cba)).append("");
		writer.append(newline);

		writer.append("seq similarity: ");
        writer.append(String.valueOf(afpChain.getSimilarity())).append("");
		writer.append(newline);

		writer.append("TM-score: ");
        writer.append(String.valueOf(afpChain.getTMScore())).append("");
		writer.append(newline);

		return writer.toString();
	}

	/**
	 * Output in FatCatCore format
	 * 
	 * <p>Note that if a circular permutation has occured the residue numbers may
	 * be innaccurate.
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @param printLegend
	 * @param longHeader
	 * @param showHTML
	 * @param showAlignmentBlock
	 * @return
	 */
	public static String toFatCatCore(
			AFPChain afpChain, 
			Atom[] ca1, 
			Atom[] ca2, 
			boolean printLegend, boolean longHeader, boolean showHTML, boolean showAlignmentBlock){

		//TODO The sequence numbers are inaccurate if a !afpChain.isSequential()

		String name1 = afpChain.getName1();
		String name2 = afpChain.getName2();
		int ca1Length = afpChain.getCa1Length();
		int ca2Length = afpChain.getCa2Length();

		int blockNum = afpChain.getBlockNum();
		int totalLenIni = afpChain.getTotalLenIni();
		double totalRmsdIni = afpChain.getTotalRmsdIni();
		int optLength = afpChain.getOptLength();
		double totalRmsdOpt = afpChain.getTotalRmsdOpt();
		double chainRmsd = afpChain.getChainRmsd();
		double alignScore = afpChain.getAlignScore();
		int alnLength = afpChain.getAlnLength();
		int gapLen = afpChain.getGapLen();
		List<AFP> afpSet = afpChain.getAfpSet();

		double similarity = afpChain.getSimilarity();
		double identity = afpChain.getIdentity();

		if (similarity <0  || identity < 0){
			afpChain.calcSimilarity();
			similarity = afpChain.getSimilarity();
			identity = afpChain.getIdentity();
		}



		String algorithmName = afpChain.getAlgorithmName();
		//String version = afpChain.getVersion();

		double probability = afpChain.getProbability();


		int afpNum = afpSet.size();

		int[] blockGap = afpChain.getBlockGap();


		double[] blockScore = afpChain.getBlockScore();
		double[] blockRmsd = afpChain.getBlockRmsd();
		int[] blockSize = afpChain.getBlockSize();


		int alnbeg1 = afpChain.getAlnbeg1();
		int alnbeg2 = afpChain.getAlnbeg2();

		char[] alnseq1 = afpChain.getAlnseq1();
		char[] alnseq2 = afpChain.getAlnseq2();
		char[] alnsymb = afpChain.getAlnsymb();

		// == end of extractation of data values from afpChain 
		////////////////////////////////

		StringBuffer txt = new StringBuffer();

		if ( longHeader) {
			txt.append(String.format("Align %s.pdb %d with %s.pdb %d", name1, ca1Length, name2, ca2Length));
		}
		else {
			txt.append(String.format("Align %s.pdb Length1: %d with %s.pdb Length2: %d", name1, ca1Length, name2, ca2Length));
		}
		txt.append(newline);
		if ( afpChain.isShortAlign()){
			txt.append("Short match");
			return txt.toString();
		}
		//txt.append(String.format("raw-score: %.2f norm.-score: %.2f ", alignScore, normAlignScore));

		if ( longHeader ) {
			txt.append(String.format( "Twists %d ini-len %d ini-rmsd %.2f opt-equ %d opt-rmsd %.2f chain-rmsd %.2f Score %.2f align-len %d gaps %d (%.2f%%)",
					blockNum - 1, totalLenIni, totalRmsdIni, optLength, totalRmsdOpt, chainRmsd, alignScore, 
					alnLength, gapLen, (100.0 * (double)gapLen/(double)alnLength)) );
			txt.append(newline);

		}  else {

			if ( ! longHeader)
				printScore(txt,algorithmName,probability,longHeader);

			printScoresInLines(afpChain, blockNum, optLength, totalRmsdOpt, alignScore, alnLength, gapLen, identity, similarity,txt);
		}


		//txt.append(String.format("P-value %.2e Afp-num %d Identity %.2f%% Similarity %.2f%% norm.-score: %.2f"+newline, probability, afpNum, identity * 100, similarity * 100, normAlignScore));

		if ( longHeader) {
			printScore(txt,algorithmName,probability,longHeader);

			txt.append(String.format("Afp-num %d Identity %.2f%% Similarity %.2f%%", afpNum, identity * 100, similarity * 100));
			txt.append(newline);
		} 

		int i;
		double gap;

		if ( longHeader ){
			int fragLen = 8 ; // FatCatParameters.DEFAULT_FRAGLEN;
			for(i = 0; i < blockNum; i ++)  {
				gap = (double)blockGap[i] /( (double)blockGap[i] + fragLen * blockSize[i]);
				txt.append(String.format( "Block %2d afp %2d score %5.2f rmsd %5.2f gap %d (%.2f%%)",
						i, blockSize[i], blockScore[i], blockRmsd[i], blockGap[i], gap));
				txt.append(newline);
			}
		}

		int     linelen = 70;
		String a;
		String b;
		String c;


		int     t = 0;
		int     ap = alnbeg1;
		int     bp = alnbeg2;
		int     k, len;

		//System.out.println(alnseq1.length + " " + alnseq1.toString());
		
		while((alnLength - t) > 0)      {
			if(alnLength - t > linelen)     len = linelen;
			else    len = alnLength - t;

			if ( ap >= ca1.length)
				break;
			if ( bp >= ca2.length)
				break;

			String pdb1 = ca1[ap].getGroup().getResidueNumber().toString();
			String pdb2 = ca2[bp].getGroup().getResidueNumber().toString();
			

			//System.err.println("t,len:"+t+":"+len);
			String lseq1 = new String(alnseq1).substring(t,t+len); 
			String lseq2 = new String(alnseq2).substring(t,t+len); 
			String lsymb = new String(alnsymb).substring(t,t+len); 

			//System.err.println("B:" + b);


			// check conservation and color accordingly, if requested by user.
			if ( showHTML ) {
				a = "";
				b = "";
				c = "";

				//	<span class=\"m\">|</span> ... Structurally equivalent and identical residues 
				//  <span class=\"sm\">:</span> ... Structurally equivalent and similar residues  
				//  <span class=\"qg\">.</span> ... Structurally equivalent, but not similar residues. 

				for (int pos = 0 ; pos < lseq1.length() ; pos ++){
					char c1 = lseq1.charAt(pos);
					char c2 = lseq2.charAt(pos);
					char cl = lsymb.charAt(pos);
					int block = -1 ;
					if ( cl != ' ') {
						try { 
							block = Integer.parseInt(cl+"");
						} catch (Exception e){
							//
						}
					}
					if ( cl != ' ' ){
						
						if ( showAlignmentBlock && block > -1 ) {
							a += "<span class=\"alignmentBlock1"+block+"\">" + c1 + "</span>";
							b += "<span class=\"alignmentBlock2"+block+"\">" + c2 + "</span>";
							c += "<span class=\"m\">" + cl + "</span>";
						} else {
							a += getPrefix(c1, c2, 0, block, false).toString() + c1 + "</span>";
							b += getPrefix(c1, c2, 1, block, false).toString() + c2 + "</span>";
							c += "<span class=\"m\">" + cl + "</span>";
						}

					} else if ( c1 != '-' && c2 != '-') {

						a += "<span class=\"sm\">" + c1 + "</span>";
						b += "<span class=\"sm\">" + c2 + "</span>";
						c += "<span class=\"sm\">" + cl + "</span>";


					} else {

						a += "<span class=\"qg\">" + c1 + "</span>";
						b += "<span class=\"qg\">" + c2 + "</span>";
						c += "<span class=\"qg\">" + cl + "</span>";

					}

					if(c1 != '-') ap ++;
					if(c2 != '-') bp ++;
				}


			} else {

				a = lseq1;
				b = lseq2;
				c = lsymb;
			}

			txt.append(newline);
			if ( longHeader )
				txt.append(String.format("%14s", " "));
			else 
				txt.append(String.format("%14s", " "));

			if (  longHeader ) {
				for(k = 10; k <= len; k += 10)
					txt.append("    .    :");
				if(k <= len + 5) txt.append("    .");

			} else {

				for(k = 10; k <= len; k += 10)
					txt.append("----+----|");
				if(k <= len + 5) txt.append("----+");


			}

			

			txt.append(newline);
			txt.append(String.format("Chain 1:%5s %s"+newline +"%14s%s"+newline+"Chain 2:%5s %s",
					pdb1, a, " ", c, pdb2, b));

			txt.append(newline);

			if ( ! showHTML){
				for(k = 0; k < len; k ++)       {
					if(a.charAt(k) != '-') ap ++;
					if(b.charAt(k) != '-') bp ++;
				}
			}
			t += len;



		}
		txt.append(newline);
		if ( printLegend ){
			if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || 
					algorithmName.equalsIgnoreCase(SmithWaterman3Daligner.algorithmName)){
				txt.append("Note: positions are from PDB; | means alignment of identical amino acids, : of similar amino acids ");

			} else {
				txt.append("Note: positions are from PDB; the numbers between alignments are block index");
			}
			txt.append(newline);

		}
		return txt.toString();

	}

	private static void printScoresInLines(AFPChain afpChain, int blockNum, int optLength, double totalRmsdOpt, double alignScore,
			int alnLength, int gapLen, double identity, double similarity, StringBuffer txt)
	{
		if ( blockNum - 1 > 0) {
			txt.append(String.format( "Twists %d ", blockNum -1 ));
			txt.append(newline);
		}

		txt.append(String.format("Equ: %d ", optLength));
		txt.append(newline);
		txt.append(String.format("RMSD: %.2f ", totalRmsdOpt));
		txt.append(newline);
		txt.append(String.format("Score: %.2f ", alignScore));
		txt.append(newline);
		txt.append(String.format("Align-len: %d ", alnLength));
		txt.append(newline);
		txt.append(String.format("Gaps: %d (%.2f%%)",
				gapLen, (100.0 * (double)gapLen/(double)alnLength)) );
		txt.append(newline);
		if ( afpChain.getTMScore() >= 0) {
			txt.append(String.format("TM-score: %.2f",afpChain.getTMScore()));
			txt.append(newline);
		}
		txt.append(newline);
		txt.append(String.format("Identity: %.2f%% ", identity * 100 ));
		txt.append(newline);
		txt.append(String.format("Similarity: %.2f%%", similarity * 100));
		txt.append(newline);
	}

	private static void printScore(StringBuffer txt,
			String algorithmName,
			double probability, 
			boolean longHeader)
	{
		if ( algorithmName.equalsIgnoreCase(CeMain.algorithmName) || algorithmName.equalsIgnoreCase(CeSideChainMain.algorithmName) ){
			txt.append(String.format("Z-score %.2f ", probability));
			if ( ! longHeader)
				txt.append(newline);
		} else if ( algorithmName.equalsIgnoreCase(SmithWaterman3Daligner.algorithmName)) {

		} else {
			if ( longHeader ){
				txt.append(String.format("P-value %.2e ",probability));
			}  else {
				txt.append(String.format("P-value: %.2e ",probability));
				txt.append(newline);
			}
		}



	}

	/**
	 * Prints the afpChain as a nicely formatted alignment, including alignment
	 * statistics, the aligned sequences themselves, and information about the 
	 * superposition.
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return a String representation as it is used on the RCSB PDB web site for display.
	 */
	public static String toWebSiteDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2 ){
		boolean showAlignmentBlock   = true;
		return toWebSiteDisplay(afpChain, ca1, ca2, showAlignmentBlock);
	}
	
	/**
	 * Prints the afpChain as a nicely formatted alignment, including alignment
	 * statistics, the aligned sequences themselves, and information about the 
	 * superposition.
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * 
	 * @return a String representation as it is used on the RCSB PDB web site for display.
	 */
	public static String toWebSiteDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2, boolean showAlignmentBlock){
		
		boolean printLegend = true;
		boolean longHeader  = true;
		boolean showHTML = true;
		
		if ( afpChain.getAlgorithmName().equalsIgnoreCase(FatCatFlexible.algorithmName)) {

		
			String msg =  toFatCatCore(afpChain,ca1,ca2,printLegend,longHeader,showHTML, showAlignmentBlock);

			return msg;
		}

		boolean showSeq = true;

		AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2, showSeq);


		//      String msg= toFatCatCore(afpChain,ca1,ca2, printLegend,longHeader);
		//


		String msg = toPrettyAlignment(afpChain, ca1, ca2, showHTML, showAlignmentBlock);


		msg = msg + newline + 
		"     <span class=\"m\">|</span> ... Structurally equivalent and identical residues " + newline +
		"     <span class=\"sm\">:</span> ... Structurally equivalent and similar residues " + newline + 
		"     <span class=\"qg\">.</span> ... Structurally equivalent, but not similar residues. " + newline;

		msg += newline;
		msg += "     To calculate the coordinates of chain 2 aligned on chain 1 apply the following transformation: ";
		msg += newline;
		msg += newline;
		msg += toRotMat(afpChain);
		return msg;


	}



	private static String toPrettyAlignment(AFPChain afpChain, Atom[] ca1, Atom[] ca2, boolean showHTML, boolean showAlignmentBlock) {
		String name1 = afpChain.getName1();
		String name2 = afpChain.getName2();
		int ca1Length = afpChain.getCa1Length();
		int ca2Length = afpChain.getCa2Length();

		int blockNum = afpChain.getBlockNum();


		int optLength = afpChain.getOptLength();
		double totalRmsdOpt = afpChain.getTotalRmsdOpt();

		double alignScore = afpChain.getAlignScore();
		int alnLength = afpChain.getAlnLength();
		int gapLen = afpChain.getGapLen();


		double similarity = afpChain.getSimilarity();
		double identity = afpChain.getIdentity();

		if (similarity <0  || identity < 0){
			afpChain.calcSimilarity();
			similarity = afpChain.getSimilarity();
			identity = afpChain.getIdentity();
		}


		String algorithmName = afpChain.getAlgorithmName();
		//String version = afpChain.getVersion();

		double probability = afpChain.getProbability();


		// == end of extractation of data values from afpChain 

		StringBuffer txt = new StringBuffer();

		txt.append(String.format("Align %s.pdb Length1: %d with %s.pdb Length2: %d", name1, ca1Length, name2, ca2Length));

		txt.append(newline);

		if ( afpChain.isShortAlign()){
			txt.append("Short match");
			return txt.toString();
		}

		printScore(txt, algorithmName, probability, false);
		printScoresInLines(afpChain, blockNum, optLength, totalRmsdOpt, alignScore, alnLength, gapLen,identity, similarity, txt);
		txt.append(newline);

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();


		int i, j,p1, p2;

		int k;
		int p1b = 0;
		int p2b = 0;

		int     len = 0;
		StringWriter alnseq1 = new StringWriter();
		StringWriter alnseq2 = new StringWriter();
		StringWriter alnsymb = new StringWriter();
		StringWriter header1  = new StringWriter();
		StringWriter footer1  = new StringWriter();
		StringWriter header2  = new StringWriter();
		StringWriter footer2  = new StringWriter();
		StringWriter block    = new StringWriter();

		int aligPos = -1;
		for(i = 0; i < blockNum; i ++)  {   

			for(j = 0; j < optLen[i]; j ++) {

				p1 = optAln[i][0][j];
				p2 = optAln[i][1][j];
				aligPos++;

				//               System.out.println(p1 + " " + p2 + " " +  footer2.toString());

				if ( len == 0){
					//the first position of sequence in alignment
					formatStartingText(p1,p2,header1,header2,footer1,footer2,ca1,ca2);
				} else {
					// check for gapped region
					int lmax = (p1 - p1b - 1)>(p2 - p2b - 1)?(p1 - p1b - 1):(p2 - p2b - 1);
					for(k = 0; k < lmax; k ++)      {


						formatGappedRegion(ca1, ca2, txt, p1, p2, k, p1b, p2b, alnseq1, alnseq2, alnsymb, header1, footer1, header2,
								footer2, block,len, showHTML);           
						len++;
						doLenCheck(len,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block, showHTML)  ;              
					}
				}

				// ALIGNED REGION
				//           System.out.println(len + " >" + header1.toString() +"< ");
				//           System.out.println(len + " >" + header2.toString() +"< ");   
				//           System.out.println(len + " >" + alnseq1.toString() +"< ");
				//           System.out.println(len + " >" + alnsymb.toString() +"< ");
				//           System.out.println(len + " >" + alnseq2.toString() +"< ");
				//           System.out.println(len + " >" + footer1.toString() +"< ");
				formatAlignedRegion(afpChain, ca1, ca2, p1, p2, alnseq1, alnseq2, alnsymb, header1, footer1, header2, footer2, block,len, aligPos, showHTML, showAlignmentBlock);
				//            System.out.println(len + " >" + header1.toString() +"< ");
				//            System.out.println(len + " >" + header2.toString() +"< ");   
				//            System.out.println(len + " >" + alnseq1.toString() +"< "); 
				//            System.out.println(len + " >" + alnsymb.toString() +"< ");
				//            System.out.println(len + " >" + alnseq2.toString() +"< ");
				//            System.out.println(len + " >" + footer1.toString() +"< ");

				len++;

				doLenCheck(len,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block, showHTML)  ;

				p1b = p1;
				p2b = p2;

				//header1.append(newline);
				//header2.append(newline);

			}

		}

		alnLength = len;

		doLenCheck(LINELENGTH,txt,header1,header2,alnseq1,alnsymb,alnseq2,footer1, footer2,block, showHTML);
		return txt.toString();
	}

	/**
	 * Prints the alignment in the simplest form: a list of aligned residues.
	 * Format is one line per residue pair, tab delimited:
	 * <ul><li>1. PDB number. Includes insertion code</li>
	 * <li>1. Chain.</li>
	 * <li>1. Amino Acid. Three letter code.</li>
	 * <li>2. PDB number.</li>
	 * <li>2. Chain.</li>
	 * <li>2. Amino Acid.</li>
	 * </ul>
	 * example:
	 * <code>152	A	ALA	161S	A	VAL</code>
	 * <p>Note that this format loses information about blocks.
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return a String representation of the aligned pairs.
	 */
	public static String toAlignedPairs(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		StringWriter pairs = new StringWriter();

		//Write structure names & PDB codes
		pairs.append("#Struct1:\t");
		pairs.append(afpChain.getName1());
		pairs.append("\n");
		pairs.append("#Struct2:\t");
		pairs.append(afpChain.getName2());
		pairs.append("\n");

		//Write optimally aligned pairs
		pairs.append("#Num1\tChain1\tAA1\tNum2\tChain2\tAA2\n");
		int[][][] optAln = afpChain.getOptAln();
		int[] blockLen = afpChain.getOptLen();
		for( int block=0;block<afpChain.getBlockNum(); block++) {
			for(int i=0;i<blockLen[block];i++) {
				Atom atom1 = ca1[ optAln[block][0][i] ];
				Atom atom2 = ca2[ optAln[block][1][i] ];

				pairs.append(atom1.getGroup().getResidueNumber().toString());
				pairs.append('\t');
				pairs.append(atom1.getGroup().getChain().getChainID());
				pairs.append('\t');
				pairs.append(atom1.getGroup().getPDBName());
				pairs.append('\t');
				pairs.append(atom2.getGroup().getResidueNumber().toString());
				pairs.append('\t');
				pairs.append(atom2.getGroup().getChain().getChainID());
				pairs.append('\t');
				pairs.append(atom2.getGroup().getPDBName());
				pairs.append('\n');
			}
		}

		return pairs.toString();
	}

	private static void formatGappedRegion(Atom[] ca1, Atom[] ca2, StringBuffer txt, int p1, int p2, int k, int p1b, int p2b,
			StringWriter alnseq1, StringWriter alnseq2, StringWriter alnsymb, StringWriter header1, StringWriter footer1,
			StringWriter header2, StringWriter footer2, StringWriter block, int len, boolean formatHTML)	{

		// DEAL WITH GAPS
		int tmppos = (p1 - p1b - 1);
		block.append("g");

		int  pos1=p1b+1+k ;
		char oneletter1 = ' ';
		try {
			oneletter1 = getOneLetter(ca1[pos1].getGroup());
		} catch (Exception e){}
		int pos2=p2b+1+k;
		char oneletter2 = ' ';
		try {
			oneletter2 = getOneLetter(ca2[pos2].getGroup());
		} catch (Exception e){}


		if(k >= tmppos) {
			//alnseq1[len] = '-';
			if (  formatHTML){
				alnseq1.append("<span class=\"qg\">-</span>");				
				header1.append(" ");
				header2.append(" ");


			} else {
				alnseq1.append('-');  
				header1.append(" ");
				header2.append(" ");
			}

		}
		else {
			if ( formatHTML){
				alnseq1.append(getPrefix(oneletter1,oneletter2,0,-1, false));
			}
			alnseq1.append(oneletter1);
			if (formatHTML){
				alnseq1.append("</span>");
			}
			formatPosition(pos1,ca1, len, header1, header2);

		}

		if(k >= (p2 - p2b - 1)) {
			//alnseq2[len] = '-';
			if ( formatHTML){
				alnseq2.append("<span class=\"qg\">-</span>");
				footer1.append(" ");
				footer2.append(" ");
			} else {
				alnseq2.append('-');
				footer1.append(" ");
				footer2.append(" ");
			}

		}
		else  {
			if ( formatHTML){
				alnseq2.append(getPrefix(oneletter1,oneletter2,1, -1, false));
			}
			alnseq2.append(oneletter2);			
			if (formatHTML){
				alnseq2.append("</span>");
			}
			formatPosition(pos2, ca2, len, footer1, footer2);

		}
		//alnsymb[len ++] = ' ';
		alnsymb.append(' ');

	}



	private static CharSequence getPrefix(char oneletter1, char oneletter2,
			int i, int blockNr, boolean showAlignmentBlock) {

		if ( oneletter1 == '-' || oneletter2 == '-' ) {
			// a gap in the alignment. 
			// label as mismatch
			return "<span class=\"qg\">";
		}

		// an aligned position
		
		if ( showAlignmentBlock && blockNr > -1){
			return "<span class=\"alignmentBlock"+(i+1)+(blockNr+1)+"\">";
		}
		
		// we return the "default" sequence alignment view...
		
		if ( oneletter1 == oneletter2)
			return "<span class=\"m\">";

		double score = AFPAlignmentDisplay.aaScore(oneletter1,oneletter2);
		if ( score > 0 )
			return "<span class=\"sm\">";

		// not similar
		return "<span class=\"qg\">";
	}

	private static void formatPosition(int pos1, Atom[] ca, int len, StringWriter header1, StringWriter header2)
	{
		int linePos = len % LINELENGTH;

		if ( header1.getBuffer().length() < linePos) {
			// fill up the buffer, we are probably shortly after the start...
			for ( int i = header1.getBuffer().length() ; i< linePos ; i++){
				header1.append(" ");
			}
		}



		Atom a = ca[pos1];
		Group g = a.getGroup();

		ResidueNumber residueNumber = g.getResidueNumber();
		pos1 = residueNumber.getSeqNum();
		boolean hasInsertionCode = false;		 
		if ( residueNumber.getInsCode() != null) {
			hasInsertionCode = true;
		} 

		if ( (pos1 %10  == 0) && ( ! hasInsertionCode)) {
			CharSequence display = getPDBPos(a);

			boolean ignoreH1 = false; 

			// make sure we don't have a problem with the left boundary...
			if ( header1.getBuffer().length()-1 > linePos) {
				ignoreH1 = true;
				System.out.println("Ignore h1: " + len + " " + header1.getBuffer().length() + " linePos: " + linePos +"  >" + header1.toString() +"<");
			}
			//System.out.println(len + " p1:" + tmp + " = " + pos1 + " " + " " + display + " " + ignoreH1);
			if ( ! ignoreH1) {
				header1.append(String.format("%-10s",display ));
				header2.append("|");
			} else {
				header2.append("|");
			}

		} else if ( hasInsertionCode){
			Character insCode = g.getResidueNumber().getInsCode();
			if ( insCode != null)
				header2.append(insCode);
			else {
				header2.append("!");
			}
		} else if ( ((pos1) %5 ) == 0 && len > 5) {
			header2.append(".");
		} else {
			if ( len > 0)
				header2.append(" ");
		}

	}



	private static void formatAlignedRegion(AFPChain afpChain, Atom[] ca1, Atom[] ca2, int p1, int p2, 
			StringWriter alnseq1, StringWriter alnseq2,
			StringWriter alnsymb, StringWriter header1, StringWriter footer1, StringWriter header2, 
			StringWriter footer2, StringWriter block, int len, int aligPos,
			boolean showHTML, boolean showAlignmentBlock)
	{
		char c1;
		char c2;
		if (( p1 < ca1.length) && (p2< ca2.length)){
			c1=  getOneLetter(ca1[p1].getGroup());
			c2 =  getOneLetter(ca2[p2].getGroup());
		} else {
			c1 = 'X';
			c2 = 'X';
		}
		
		int blockPos = -1;
		if ( afpChain.getBlockNum() > 0) {
			blockPos = AFPAlignmentDisplay.getBlockNrForAlignPos(afpChain, aligPos);
		}	
		
		
		
		double score = AFPAlignmentDisplay.aaScore(c1,c2);

		if ( showHTML) {
			alnseq1.append(getPrefix(c1,c2,  0, blockPos, showAlignmentBlock));
			alnseq2.append(getPrefix(c1,c2,  1, blockPos, showAlignmentBlock));
		}

		alnseq1.append(c1);              
		alnseq2.append(c2);

		if ( showHTML){
			alnseq1.append("</span>");
			alnseq2.append("</span>");
		}

		if ( c1 == c2){
			if ( showHTML){
				
				alnsymb.append("<span class=\"m\">|</span>");
				
			} else {
				alnsymb.append('|');
			}
			//alnsymb[len ++] = '|';
		} else {


			if ( score > 1) {
				if ( showHTML){
					alnsymb.append( "<span class=\"sm\">:</span>");
				} else {
					alnsymb.append( ':');
				}
			}
			else {
				if ( showHTML)
					alnsymb.append( "<span class=\"qg\">.</span>");					
				else
					alnsymb.append( '.');
			}
		}

		if ( p1 < ca1.length)
			formatPosition(p1, ca1,len, header1, header2);

		if ( p2 < ca2.length)
			formatPosition(p2,ca2,len, footer1, footer2);

	}

	private static void formatStartingText(int p1, int p2, StringWriter header1, StringWriter header2, StringWriter footer1,
			StringWriter footer2, Atom[] ca1, Atom[] ca2)
	{

		header1.append(String.format("%-10s", getPDBPos(ca1[p1])));
		header2.append("|");
		footer1.append(String.format("%-10s", getPDBPos(ca2[p2])));
		footer2.append("|");


	}

	private static boolean doLenCheck(int len, StringBuffer txt, StringWriter header1, StringWriter header2, StringWriter alnseq1,
			StringWriter alnsymb, StringWriter alnseq2, StringWriter footer1, StringWriter footer2, StringWriter block, boolean formatHTML)
	{

		if ( len % LINELENGTH  == 0) {

			//txt.append("|");
			txt.append(header1);
			//txt.append("|");
			txt.append(newline);
			//txt.append("|");
			txt.append(header2);
			//txt.append("|");
			txt.append(newline);
			//txt.append("|");
			txt.append(alnseq1);
			//txt.append("|");
			txt.append(newline);

			//txt.append("|");
			txt.append(alnsymb);
			//         txt.append(newline);
			//         txt.append(block);
			//txt.append("|");
			txt.append(newline);
			//txt.append("|");
			txt.append(alnseq2);
			//txt.append("|");
			txt.append(newline);
			//txt.append("|");
			txt.append(footer2);
			//txt.append("|");
			txt.append(newline);
			//txt.append("|");
			txt.append(footer1);
			//txt.append("|");
			txt.append(newline);
			txt.append(newline);
			txt.append(newline);

			if (formatHTML ) {

				int len1 = alnseq1.getBuffer().length();
				int len2 = alnseq2.getBuffer().length();
				int lens = alnsymb.getBuffer().length();
				alnseq1.getBuffer().replace(0, len1, "");
				alnseq2.getBuffer().replace(0, len2, "");         
				alnsymb.getBuffer().replace(0, lens, "");              



				header1.getBuffer().replace(0, LINELENGTH, "");
				header2.getBuffer().replace(0, LINELENGTH , "");
				footer1.getBuffer().replace(0, LINELENGTH, "");         
				footer2.getBuffer().replace(0, LINELENGTH, "");
				block.getBuffer().replace(0, LINELENGTH, "");
			} else {
				alnseq1.getBuffer().replace(0, LINELENGTH, "");
				alnseq2.getBuffer().replace(0, LINELENGTH, "");         
				alnsymb.getBuffer().replace(0, LINELENGTH, "");              
				header1.getBuffer().replace(0, LINELENGTH, "");
				header2.getBuffer().replace(0, LINELENGTH , "");
				footer1.getBuffer().replace(0, LINELENGTH, "");         
				footer2.getBuffer().replace(0, LINELENGTH, "");
				block.getBuffer().replace(0, LINELENGTH, "");
			}
			StringBuffer buf = header1.getBuffer();
			for ( int i=0;i<buf.length();i++){
				char c = buf.charAt(i);
				if ( c != ' '){
					buf.setCharAt(i, ' ');
				}
			}
			buf = footer1.getBuffer();
			for ( int i=0;i<buf.length();i++){
				char c = buf.charAt(i);
				if ( c != ' '){
					buf.setCharAt(i, ' ');
				}
			}

			return true;
		}

		return false;


	}

	private static CharSequence getPDBPos(Atom atom)
	{

		Group g = atom.getGroup();
		if ( g!= null){
			Chain c = g.getChain();
			if (c != null){
				return g.getResidueNumber().toString()+":" + c.getChainID() ;
				//return g.getPDBCode()+":" + c.getName() + "." + getOneLetter(g) ; 
			}
		}
		return "!";
	}

	private static char getOneLetter(Group g){

		try {
			Character c = StructureTools.get1LetterCode(g.getPDBName());
			return c;
		} catch (Exception e){
			return 'X';
		}
	}


	public static String toDBSearchResult(AFPChain afpChain)
	{
		StringBuffer str = new StringBuffer();

		str.append(afpChain.getName1());
		str.append("\t");
		str.append(afpChain.getName2());
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getAlignScore()));
		str.append("\t");     
		if ( afpChain.getAlgorithmName().equalsIgnoreCase(CeMain.algorithmName)){
			str.append(String.format("%.2f",afpChain.getProbability()));
		} else {
			str.append(String.format("%.2e",afpChain.getProbability()));
		}
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getTotalRmsdOpt()));
		str.append("\t");
		str.append(afpChain.getCa1Length());
		str.append("\t");
		str.append(afpChain.getCa2Length());      
		str.append("\t");
		str.append(afpChain.getCoverage1());
		str.append("\t");
		str.append(afpChain.getCoverage2());
		str.append("\t");
		str.append(String.format("%.2f",afpChain.getIdentity()));
		str.append("\t");
		str.append(afpChain.getDescription2());
		str.append("\t");
		str.append(newline);

		return str.toString();
	}

	public static String toRotMat(AFPChain afpChain)
	{

		Matrix[] blockRotationMatrix = afpChain.getBlockRotationMatrix();
		int blockNum = afpChain.getBlockNum();
		Atom[] blockShiftVector = afpChain.getBlockShiftVector();

		StringBuffer txt = new StringBuffer();

		if ( blockRotationMatrix == null || blockRotationMatrix.length < 1)
			return "";
	

		for ( int blockNr = 0 ; blockNr < blockNum  ; blockNr++){
			Matrix m = blockRotationMatrix[blockNr];
			Atom shift   = blockShiftVector[blockNr];
			if ( blockNum > 1) {
				txt.append("Operations for block " );
				txt.append(blockNr);
				txt.append(newline);
			}

			String origString = "orig";
			if ( blockNr > 0)
				origString = (blockNr)+""; 


			txt.append(String.format("     X"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,0),m.get(1,0), m.get(2,0), shift.getX()));
			txt.append( newline); 
			txt.append(String.format("     Y"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,1),m.get(1,1), m.get(2,1), shift.getY()));
			txt.append( newline);
			txt.append(String.format("     Z"+(blockNr+1)+" = (%9.6f)*X"+ origString +" + (%9.6f)*Y"+ origString +" + (%9.6f)*Z"+ origString +" + (%12.6f)",m.get(0,2),m.get(1,2), m.get(2,2), shift.getZ()));
			txt.append(newline);
		}
		return txt.toString();
	}

	public static String toCE(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
	{



		String name1 = afpChain.getName1();
		String name2 = afpChain.getName2();

		int optLength = afpChain.getOptLength();
		double totalRmsdOpt = afpChain.getTotalRmsdOpt();

		int alnLength = afpChain.getAlnLength();
		int gapLen = afpChain.getGapLen();


		double similarity = afpChain.getSimilarity();
		double identity = afpChain.getIdentity();
		if (similarity <0 || identity <0  ){
			afpChain.calcSimilarity();
			similarity = afpChain.getSimilarity();
			identity = afpChain.getIdentity();
		}

	

		double probability = afpChain.getProbability();


		int alnbeg1 = afpChain.getAlnbeg1();
		int alnbeg2 = afpChain.getAlnbeg2();

		char[] alnseq1 = afpChain.getAlnseq1();
		char[] alnseq2 = afpChain.getAlnseq2();


		long calculationTime = afpChain.getCalculationTime();

		// == end of extractation of data values from afpChain 



		StringBuffer txt = new StringBuffer();

		txt.append("Chain 1: ");
		txt.append(name1);
		txt.append(" (Size=");
		txt.append(ca1.length);
		txt.append(")");
		txt.append(newline);
		txt.append("Chain 2: ");
		txt.append(name2);
		txt.append(" (Size=");
		txt.append(ca2.length);
		txt.append(")");
		txt.append(newline);
		txt.append(newline);
		txt.append(String.format("Alignment length = %d Rmsd = %.2fA Z-Score = %.1f",optLength,totalRmsdOpt,probability));
		txt.append(String.format(" Gaps = %d(%.1f%%) CPU = %d ms. Sequence identities = %.1f%%",gapLen,( gapLen*100.0/optLength),calculationTime,identity*100));

		int     linelen = 70;
		String a;
		String b;



		int     t = 0;
		int     ap = alnbeg1;
		int     bp = alnbeg2;
		int     k, len;

		while((alnLength - t) > 0)      {
			if(alnLength - t > linelen)     len = linelen;
			else    len = alnLength - t;


			//System.err.println("t,len:"+t+":"+len);
			a = new String(alnseq1).substring(t,t+len);
			b = new String(alnseq2).substring(t,t+len);

			//System.err.println("B:" + b);

			/*
            txt.append(newline);
            txt.append(String.format("%14s", " "));

            for(k = 10; k <= len; k += 10)
                txt.append("    .    :");
            if(k <= len + 5) txt.append("    .");
			 */

			//String pdb1 = ca1[ap].getParent().getPDBCode();
			//String pdb2 = ca2[bp].getParent().getPDBCode();
			txt.append(newline);
			txt.append(String.format("Chain 1:%5s %s"+newline+"Chain 2:%5s %s",
					(ap+1), a, (bp+1), b));
			txt.append(newline);
			for(k = 0; k < len; k ++)       {
				if(a.charAt(k) != '-') ap ++;
				if(b.charAt(k) != '-') bp ++;
			}
			t += len;

		}
		txt.append(newline);

		txt.append(toRotMat(afpChain));

		return txt.toString();


	}


	
	


}
