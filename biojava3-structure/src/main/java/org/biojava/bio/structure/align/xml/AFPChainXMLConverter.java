package org.biojava.bio.structure.align.xml;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.core.util.PrettyXMLWriter;



public class AFPChainXMLConverter {


	/** Convert an afpChain to a simple XML representation
	 * 
	 * @param afpChain
	 * @return XML representation of the AFPCHain
	 */
	public synchronized static String toXML(AFPChain afpChain, Atom[] ca1, Atom[]ca2) throws IOException{
		StringWriter result = new StringWriter();	
		toXML(afpChain,result,ca1,ca2);
		return result.toString();

	}


	/** Write the XML representation to a StringWriter
	 * 
	 * @param afpChain
	 * @param swriter
	 * @throws IOException
	 */
	public synchronized static void toXML(AFPChain afpChain, StringWriter swriter,Atom[] ca1, Atom[]ca2) throws IOException{

		PrintWriter writer = new PrintWriter(swriter);
		PrettyXMLWriter xml = new PrettyXMLWriter(writer);
		
		xml.openTag("AFPChain");
		
		printXMLHeader(xml,afpChain);


		// that is the initial alignment...
		// we don't serialize that at the present.
		//int[] blockResSize = afpChain.getBlockResSize();
		//int[][][] blockResList = afpChain.getBlockResList();


		// get the alignment blocks
		int blockNum = afpChain.getBlockNum();
		//int[] optLen       = afpChain.getOptLen();
		//int[] blockSize    = afpChain.getBlockSize();
		for(int bk = 0; bk < blockNum; bk ++) {
						
			xml.openTag("block");

			printXMLBlockHeader(xml,afpChain, bk);

			if ( ca1 == null || ca2 == null) {
				try {
					printXMLEQRKnownPositions(xml,afpChain,bk);
				} catch (StructureException ex ){
					throw new IOException(ex.getMessage());
				}
			}
			else 
				printXMLEQRInferPositions(xml, afpChain,bk,ca1,ca2);

			printXMLMatrixShift(xml, afpChain, bk);

			xml.closeTag("block");
		}

		xml.closeTag("AFPChain");

		writer.close();


	}



	private static void printXMLEQRKnownPositions(PrettyXMLWriter xml,
			AFPChain afpChain, int blockNr) throws IOException, StructureException{
		int[] optLen       = afpChain.getOptLen();


		String[][][] pdbAln = afpChain.getPdbAln();
		if ( pdbAln == null){
			throw new StructureException("Can't convert to XML without known the PDB coordinates. Please provide Ca atoms and call toXML(afpChain,ca1, ca2)");
		}

		for ( int eqrNr = 0 ; eqrNr < optLen[blockNr] ; eqrNr++ ){
			
			String pdbResnum1 = pdbAln[blockNr][0][eqrNr];
			String pdbResnum2 = pdbAln[blockNr][1][eqrNr];
			
			//System.out.println(eqrNr + " got resnum: " + pdbResnum1 + " " + pdbResnum2);
			
			String[] spl1 = pdbResnum1.split(":");
			String[] spl2 = pdbResnum2.split(":");

			String chain1 = spl1[0];
			String pdbres1 = spl1[1];

			String chain2 = spl2[0];
			String pdbres2 = spl2[1];

			xml.openTag("eqr");
			xml.attribute("eqrNr",eqrNr+"");
			xml.attribute("pdbres1",pdbres1);
			xml.attribute("chain1", chain1);
			xml.attribute("pdbres2",pdbres2);
			xml.attribute("chain2", chain2);
			xml.closeTag("eqr");
		}
	}


	public static void printXMLEQRInferPositions(PrettyXMLWriter xml,			
			AFPChain afpChain, int bk, Atom[] ca1, Atom[] ca2)  throws IOException{
		
		int[] optLen       = afpChain.getOptLen();
		
		if ( optLen == null)
			return;
				
		int[][][] optAln   = afpChain.getOptAln();
		
		
		for ( int pos=0;pos< optLen[bk];pos++){
			int pos1 = optAln[bk][0][pos];
			int pos2 = optAln[bk][1][pos];
			xml.openTag("eqr");
			xml.attribute("eqrNr",pos+"");
			xml.attribute("pdbres1",ca1[pos1].getGroup().getResidueNumber().toString());
			xml.attribute("chain1", ca1[pos1].getGroup().getChain().getChainID());
			xml.attribute("pdbres2",ca2[pos2].getGroup().getResidueNumber().toString());
			xml.attribute("chain2", ca2[pos2].getGroup().getChain().getChainID());
			
			xml.closeTag("eqr");
			//System.out.println("aligned position: " + pos1  + ":" + pos2 + 
			//" pdbresnum " + ca1[pos1].getGroup().getResidueNumber().toString() + " " +
			//ca1[pos1].getParent().getPDBName()+":" + 
			//ca2[pos2].getGroup().getResidueNumber().toString() + " " + ca2[pos2].getParent().getPDBName());
	 
		}	 

	}


	private static void printXMLBlockHeader(PrettyXMLWriter xml,
			AFPChain afpChain,int blockNr) throws IOException{

		int bk = blockNr;
		int[] blockSize    = afpChain.getBlockSize();
		//if ( blockSize[bk] == 0) {
		//	return;
		//}
		int[] blockGap     = afpChain.getBlockGap();
		
		double[]blockScore = afpChain.getBlockScore();
		double[] blockRmsd = afpChain.getBlockRmsd();

		xml.attribute("blockNr",bk+"");
		xml.attribute("blockSize", blockSize[bk]+"");
		xml.attribute("blockScore", String.format("%5.2f",blockScore[bk]).trim());
		xml.attribute("blockRmsd", String.format("%5.2f",blockRmsd[bk]).trim());
		xml.attribute("blockGap", blockGap[bk]+"");

	}


	private static void printXMLMatrixShift(PrettyXMLWriter xml,
			AFPChain afpChain, int blockNr)  throws IOException {
	   
	   Matrix[] ms     = afpChain.getBlockRotationMatrix();
	   if ( ms == null || ms.length == 0)
          return;
	   
	   Matrix matrix = ms[blockNr];	
	   if ( matrix == null)
		   return;
	   xml.openTag("matrix");
		
			
		for (int x=0;x<3;x++){
			for (int y=0;y<3;y++){
				String key = "mat"+(x+1)+(y+1);
				xml.attribute(key,String.format("%.6f",matrix.get(x,y)));
			}
		}
		xml.closeTag("matrix");

		Atom[]   shifts = afpChain.getBlockShiftVector();
		Atom shift = shifts[blockNr];
		xml.openTag("shift");
		xml.attribute("x", String.format("%.3f",shift.getX()));
		xml.attribute("y", String.format("%.3f",shift.getY()));
		xml.attribute("z", String.format("%.3f",shift.getZ()));
		xml.closeTag("shift");

	}


	public static String toXML(AFPChain afpChain) throws IOException{

		return toXML(afpChain, null, null);
	}


	public static void printXMLHeader(PrettyXMLWriter xml, AFPChain afpChain) throws IOException{
		xml.attribute("name1", afpChain.getName1());
		xml.attribute("name2", afpChain.getName2());
		xml.attribute("method", afpChain.getAlgorithmName());
		xml.attribute("version" , afpChain.getVersion());
		xml.attribute("alnLength", afpChain.getAlnLength() + "");
		xml.attribute("blockNum", afpChain.getBlockNum() + "");
		xml.attribute("gapLen", afpChain.getGapLen() + "");
		xml.attribute("optLength", afpChain.getOptLength() + "");
		xml.attribute("totalLenIni", afpChain.getTotalLenIni() + "");

		xml.attribute("alignScore", String.format("%5.2f", afpChain.getAlignScore() ).trim());
		xml.attribute("chainRmsd",  String.format("%5.2f", afpChain.getChainRmsd() ).trim());
		xml.attribute("identity",String.format("%5.4f", afpChain.getIdentity() ).trim());
		xml.attribute("normAlignScore", String.format("%5.2f",afpChain.getNormAlignScore()).trim());
		xml.attribute("probability", String.format("%.2e", afpChain.getProbability() ).trim());
		xml.attribute("similarity", String.format("%5.4f", afpChain.getSimilarity() ).trim());
		
		xml.attribute("similarity1", afpChain.getCoverage1() + "");
		xml.attribute("similarity2", afpChain.getCoverage2() + "");
		xml.attribute("totalRmsdIni", String.format("%5.2f",afpChain.getTotalRmsdIni() ).trim());
		xml.attribute("totalRmsdOpt", String.format("%5.2f",afpChain.getTotalRmsdOpt() ).trim());
		xml.attribute("ca1Length", afpChain.getCa1Length()+"");
		xml.attribute("ca2Length", afpChain.getCa2Length()+"");
		xml.attribute("afpNum",afpChain.getAfpSet().size()+"");
		xml.attribute("alignScoreUpdate",String.format("%5.2f",afpChain.getAlignScoreUpdate()).trim());
		xml.attribute("time", String.format("%d",afpChain.getCalculationTime()));
		if ( afpChain.getTMScore() != -1){
			xml.attribute("tmScore", String.format("%.2f",afpChain.getTMScore()));
		}
		
		// test if alignment is CP:
		if ( ! AlignmentTools.isSequentialAlignment(afpChain,false)) {
			xml.attribute("cp","true");
		}
	}
}