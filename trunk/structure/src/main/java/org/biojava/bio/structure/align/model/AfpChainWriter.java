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

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;

public class AfpChainWriter
{

	public static final String newline = System.getProperty("line.separator");

	public static String toFatCat(AFPChain afpChain, Atom[] ca1, Atom[] ca2)
	{

		boolean printLegend = true;
		boolean longHeader  = true;

		return toFatCatCore(afpChain, ca1, ca2, printLegend, longHeader);   
	}

	public static String toFatCatCore(
			AFPChain afpChain, 
			Atom[] ca1, 
			Atom[] ca2, 
			boolean printLegend, boolean longHeader){

		if(!afpChain.isSequentialAlignment()) {
			return "Can't display circular permutations";
		}

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
		}


		//txt.append(String.format("P-value %.2e Afp-num %d Identity %.2f%% Similarity %.2f%% norm.-score: %.2f"+newline, probability, afpNum, identity * 100, similarity * 100, normAlignScore));

		if ( longHeader)
			printScore(txt,algorithmName,probability,longHeader);



		if (  longHeader) {

			txt.append(String.format("Afp-num %d Identity %.2f%% Similarity %.2f%%", afpNum, identity * 100, similarity * 100));
			txt.append(newline);
		} else {

			txt.append(String.format("Identity: %.2f%% ", identity * 100 ));
			txt.append(newline);
			txt.append(String.format("Similarity: %.2f%%", similarity * 100));
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
		while((alnLength - t) > 0)      {
			if(alnLength - t > linelen)     len = linelen;
			else    len = alnLength - t;


			//System.err.println("t,len:"+t+":"+len);
			a = new String(alnseq1).substring(t,t+len);
			b = new String(alnseq2).substring(t,t+len);
			c = new String(alnsymb).substring(t,t+len);
			//System.err.println("B:" + b);

			txt.append(newline);
			txt.append(String.format("%14s", " "));

         if ( ! longHeader ) {
			for(k = 10; k <= len; k += 10)
				txt.append("    .    :");
			if(k <= len + 5) txt.append("    .");

         } else {
            
            for(k = 10; k <= len; k += 10)
               txt.append("----+----|");
            if(k <= len + 5) txt.append("----+");
         }
         
         
			String pdb1 = ca1[ap].getParent().getPDBCode();
			String pdb2 = ca2[bp].getParent().getPDBCode();
			txt.append(newline);
			txt.append(String.format("Chain 1:%5s %s"+newline +"%14s%s"+newline+"Chain 2:%5s %s",
					pdb1, a, " ", c, pdb2, b));
			txt.append(newline);
			for(k = 0; k < len; k ++)       {
				if(a.charAt(k) != '-') ap ++;
				if(b.charAt(k) != '-') bp ++;
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

	public static String toWebSiteDisplay(AFPChain afpChain, Atom[] ca1, Atom[] ca2){
		if ( afpChain.getAlgorithmName().equalsIgnoreCase("jFatCat_flexible")) {
			String msg =  toFatCat(afpChain,ca1,ca2) ;

			return msg;
		}

		boolean showSeq = true;

		AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2, showSeq);

		boolean printLegend = false;
		boolean longHeader  = false;

		String msg= toFatCatCore(afpChain,ca1,ca2, printLegend,longHeader);

		msg = msg + newline + 
		"     | ... Structurally equivalend and identical residues " + newline +
		"     : ... Structurally equivalend and similar residues " + newline + 
		"     . ... Structurally equivalent, but not similar residues. " + newline;

		msg += newline;
		msg += "     To calculate the coordinates of chain 2 aligned on chain 1 apply the following transformation: ";
		msg += newline;
		msg += newline;
		msg += toRotMat(afpChain);
		return msg;


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
		str.append(afpChain.getSimilarity1());
		str.append("\t");
		str.append(afpChain.getSimilarity2());
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
		if (similarity == -1 || identity == -1){
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
