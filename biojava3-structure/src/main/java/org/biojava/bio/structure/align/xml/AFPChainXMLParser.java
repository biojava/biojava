package org.biojava.bio.structure.align.xml;


import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;

import org.w3c.dom.Document;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

//http://www.developerfusion.com/code/2064/a-simple-way-to-read-an-xml-file-in-java/

public class AFPChainXMLParser
{

	
	/** new utility method that checks that the order of the pair in the XML alignment is correct and flips the direction if needed
	 * 
	 * @param xml
	 * @param name1
	 * @param name1
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	 public static AFPChain fromXML(String xml, String name1, String name2, Atom[] ca1, Atom[] ca2) throws StructureException{
			AFPChain[] afps = parseMultiXML( xml);
			if ( afps.length > 0 ) {

				AFPChain afpChain = afps[0];
				
				String n1 = afpChain.getName1();
				String n2 = afpChain.getName2();
				
				if ( n1 == null )
					n1 = "";
				if ( n2 == null)
					n2 = "";
				
				//System.out.println("from AFPCHAIN: " + n1 + " " + n2);
				if ( n1.equals(name2) && n2.equals(name1)){
					// flipped order
					//System.out.println("AfpChain in wrong order, flipping...");
					afpChain  = AFPChainFlipper.flipChain(afpChain);
				}
				rebuildAFPChain(afpChain, ca1, ca2);

				return afpChain;
			}
			return null;
		 
	 }
	
   public static AFPChain fromXML(String xml, Atom[] ca1, Atom[] ca2)
	{
		AFPChain[] afps = parseMultiXML( xml);
		if ( afps.length > 0 ) {

			AFPChain afpChain = afps[0];
			rebuildAFPChain(afpChain, ca1, ca2);

			return afpChain;
		}
		return null;
	}

   /** returns true if the alignment XML contains an error message
    * 
    * @param xml
    * @return flag if there was an Error while processing the alignment. 
    */
   public static boolean isErrorXML(String xml){
      
      if ( xml.contains("error=\""))
         return true;
      
      return false;
         
      
   }
   
	/** Takes an XML representation of the alignment and flips the positions of name1 and name2
	 * 
	 * @param xml String representing the alignment
	 * @return XML representation of the flipped alignment
	 */
	public static String flipAlignment(String xml) throws IOException,StructureException{
		AFPChain[] afps = parseMultiXML( xml);
		if ( afps.length < 1 ) 
			return null;

		if ( afps.length == 1) {
			AFPChain newChain = AFPChainFlipper.flipChain(afps[0]);
			if ( newChain.getAlgorithmName() == null) {
			   newChain.setAlgorithmName(AFPChain.DEFAULT_ALGORITHM_NAME);
			}
			return AFPChainXMLConverter.toXML(newChain);
		}
		throw new StructureException("not Implemented yet!");
	}


	/**  replace the PDB res nums with atom positions:
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 */
	public static void rebuildAFPChain(AFPChain afpChain, Atom[] ca1, Atom[] ca2){

	   if ( afpChain.getAlgorithmName() == null) {
	      afpChain.setAlgorithmName(AFPChain.DEFAULT_ALGORITHM_NAME);	      
	   }
	   if ( afpChain.getVersion() == null){
	      afpChain.setVersion("1.0");
	   }
	   
		int blockNum  = afpChain.getBlockNum();
		int ca1Length = afpChain.getCa1Length();
		int ca2Length = afpChain.getCa2Length();

		int minLength = Math.min(ca1Length, ca2Length);
		int[][][] optAln = new int[blockNum][2][minLength];
				
		int[][][] blockResList = afpChain.getBlockResList();
		if ( blockResList == null){
		   blockResList = new int[blockNum][2][minLength];
		}
		int[] optLen = afpChain.getOptLen();

		String[][][] pdbAln = afpChain.getPdbAln();
		int[] verifiedOptLen = null;
		if ( optLen != null)
		  verifiedOptLen = afpChain.getOptLen().clone();
		else {
			System.err.println("did not find optimal alignment, building up empty alignment.");
			optLen = new int[1];
			optLen[0] = 0;
		}
		for (int blockNr = 0 ; blockNr < blockNum ; blockNr++){

			//System.out.println("got block " + blockNr + " size: " + optLen[blockNr]);
		   int verifiedEQR = -1;
			for ( int eqrNr = 0 ; eqrNr < optLen[blockNr] ; eqrNr++ ){
				String pdbResnum1 = pdbAln[blockNr][0][eqrNr];
				String pdbResnum2 = pdbAln[blockNr][1][eqrNr];

				//System.out.println(blockNr + " " + eqrNr + " got resnum: " + pdbResnum1 + " " + pdbResnum2);
				String[] spl1 = pdbResnum1.split(":");
				String[] spl2 = pdbResnum2.split(":");

				String chain1 = spl1[0];
				String pdbres1 = spl1[1];

				String chain2 = spl2[0];
				String pdbres2 = spl2[1];

				int pos1 = getPositionForPDBresunm(pdbres1,chain1,ca1);
				int pos2 = getPositionForPDBresunm(pdbres2,chain2,ca2);
				
				if ( pos1 == -1 || pos2 == -1 ){
				   // this can happen when parsing old files that contained Calcium atoms...
				   System.err.println("AFPChainXMLParser: warning: pos1: " +pos1 + " " + pdbResnum1 + " pos2: " + pos2 + " " + pdbResnum2 +  " should never be -1. Probably parsing an.");
				   verifiedOptLen[blockNr]-- ;
				   continue;
				}

				verifiedEQR++;
				//System.out.println(blockNr + " " + eqrNr + " " + pos1 + " " + pos2);
				optAln[blockNr][0][verifiedEQR] = pos1;
				optAln[blockNr][1][verifiedEQR] = pos2;
				blockResList[blockNr][0][verifiedEQR] = pos1;
				blockResList[blockNr][1][verifiedEQR] = pos2;
			}
		}

		afpChain.setOptLen(verifiedOptLen);
		afpChain.setOptAln(optAln);
		afpChain.setBlockResList(blockResList);
		// build up alignment image:
		AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2);


	}

	public static AFPChain[] parseMultiXML(String xml) {
		List<AFPChain> afpChains = new ArrayList<AFPChain>();

		try
		{
			//Convert string to XML document
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = factory.newDocumentBuilder();
			InputSource inStream = new InputSource();
			inStream.setCharacterStream(new StringReader(xml));
			Document doc = db.parse(inStream);

			// normalize text representation
			doc.getDocumentElement().normalize();


			//Element rootElement = doc.getDocumentElement();

			NodeList listOfAFPChains = doc.getElementsByTagName("AFPChain");
			//int numArrays = listOfArrays.getLength();			
			// go over the blocks
			for(int afpPos=0; afpPos<listOfAFPChains.getLength() ; afpPos++)
			{

				AFPChain a = new AFPChain();
				a.setAlgorithmName(AFPChain.DEFAULT_ALGORITHM_NAME);
				a.setVersion("1.0");
				Node rootElement       = listOfAFPChains.item(afpPos);

				a.setName1(getAttribute(rootElement,"name1"));				
				a.setName2(getAttribute(rootElement,"name2"));
				String algoname = getAttribute(rootElement,"method");
				if ( algoname != null) {
					a.setAlgorithmName(algoname);					
				}
				String version = getAttribute(rootElement,"version");
				if ( version != null)
					a.setVersion(version);

				a.setAlnLength(	new Integer(getAttribute(rootElement,"alnLength")).intValue());
				a.setBlockNum(		new Integer(getAttribute(rootElement,"blockNum")).intValue());
				a.setGapLen(		new Integer(getAttribute(rootElement,"gapLen")).intValue());
				a.setOptLength(	new Integer(getAttribute(rootElement,"optLength")).intValue());
				a.setTotalLenIni(	new Integer(getAttribute(rootElement,"totalLenIni")).intValue());
				a.setBlockNum(		new Integer(getAttribute(rootElement,"blockNum")).intValue());

				if ( a.getAlgorithmName().equals(CeCPMain.algorithmName)){
                   a.setSequentialAlignment(a.getBlockNum() == 1);
                }
				
				a.setAlignScore(new Double(getAttribute(rootElement,"alignScore")).doubleValue());
				a.setChainRmsd(new Double(getAttribute(rootElement,"chainRmsd")).doubleValue());
				Double identity = new Double(getAttribute(rootElement,"identity")).doubleValue();				
				a.setIdentity(identity);
				
				a.setNormAlignScore(new Double(getAttribute(rootElement,"normAlignScore")).doubleValue());
				a.setProbability(new Double(getAttribute(rootElement,"probability")).doubleValue());
				a.setSimilarity(new Double(getAttribute(rootElement,"similarity")).doubleValue());
				a.setTotalRmsdIni(new Double(getAttribute(rootElement,"totalRmsdIni")).doubleValue());
				a.setTotalRmsdOpt(new Double(getAttribute(rootElement,"totalRmsdOpt")).doubleValue());
				a.setAlignScoreUpdate(new Double(getAttribute(rootElement,"alignScoreUpdate")).doubleValue());
				int ca1Length = new Integer(getAttribute(rootElement,"ca1Length")).intValue();
				a.setCa1Length(ca1Length);
				int ca2Length = new Integer(getAttribute(rootElement,"ca2Length")).intValue();						
				a.setCa2Length(ca2Length);

				String tmScoreS = getAttribute(rootElement,"tmScore");
				if ( tmScoreS != null) {
					Double tmScore = null;
					try {
					 tmScore = Double.parseDouble(tmScoreS);
					} catch (Exception e){						
					}
					a.setTMScore(tmScore);
				}
				
				String calcTimeS = getAttribute(rootElement,"time");
				Long calcTime = -1L;
				if ( calcTimeS != null){
				
					try {
						calcTime = Long.parseLong(calcTimeS);
						
					} catch (Exception e){
						e.printStackTrace();
					}					
				}
				a.setCalculationTime(calcTime);
				
				Matrix[] ms = new Matrix[a.getBlockNum()];
				a.setBlockRotationMatrix(ms);
				Atom[] blockShiftVector = new Atom[a.getBlockNum()];
				a.setBlockShiftVector(blockShiftVector);

				int afpNum = new Integer(getAttribute(rootElement,"afpNum")).intValue();
				List<AFP> afpSet = new ArrayList<AFP>();
				for (int afp=0;afp<afpNum;afp++){
					afpSet.add( new AFP());	
				}

				a.setAfpSet(afpSet);

				int minLength = Math.min(ca1Length, ca2Length);
				a.setFocusRes1(new int[minLength]);
				a.setFocusRes2(new int[minLength]);


				//NodeList listOfBlocks = doc.getElementsByTagName("block");
				NodeList listOfBlocks = rootElement.getChildNodes();

				//int numArrays = listOfArrays.getLength();

				// go over the blocks
				for(int i=0; i<listOfBlocks.getLength() ; i++)
				{
					Node block       = listOfBlocks.item(i);

					// we only look at blocks.
					if (! block.getNodeName().equals("block"))
						continue;

					processBlock(block, a, minLength);


				}

				afpChains.add(a);
			}
		}
		catch (SAXParseException err) 
		{
			System.out.println ("** Parsing error" + ", line " 
					+ err.getLineNumber () + ", uri " + err.getSystemId ());
			System.out.println(" " + err.getMessage ());
		}
		catch (SAXException e)
		{
			Exception x = e.getException ();
			((x == null) ? e : x).printStackTrace ();
		}
		catch (Throwable t)
		{
			t.printStackTrace ();
		}

		return afpChains.toArray(new AFPChain[afpChains.size()]);
	}


	private static  void processBlock(Node block, AFPChain a, int minLength){
		NodeList valList = block.getChildNodes();
		int numChildren  = valList.getLength();

		NamedNodeMap map = block.getAttributes();

		int blockNum = a.getBlockNum();

		int[] 	  optLen 			= a.getOptLen();
		if ( optLen == null ) 
			optLen = new int[blockNum];

		String[][][] pdbAln = a.getPdbAln();
		if ( pdbAln == null)
			pdbAln         = new String[blockNum][2][minLength];

		//int[][][] optAln 			= new int[blockNum][2][minLength];
		int[]     blockGap = a.getBlockGap();
		if ( blockGap == null )
			blockGap = new int[blockNum];
		int[]     blockSize= a.getBlockSize();
		if ( blockSize == null)
			blockSize = new int[blockNum];

		double[]  blockScore = a.getBlockScore();
		if ( blockScore == null)
			blockScore = new double[blockNum];
		double[]  blockRmsd = a.getBlockRmsd();
		if (blockRmsd == null )
			blockRmsd = new double[blockNum];
		Matrix[] ms     = a.getBlockRotationMatrix();
		Atom[] shifts = a.getBlockShiftVector();

		int blockNr = new Integer( map.getNamedItem("blockNr").getTextContent()).intValue();

		int thisBlockGap = new Integer(map.getNamedItem("blockGap").getTextContent()).intValue();
		blockGap[blockNr] = thisBlockGap;

		int thisBlockSize = new Integer(map.getNamedItem("blockSize").getTextContent()).intValue();
		blockSize[blockNr] = thisBlockSize;

		double thisBlockScore = new Double(map.getNamedItem("blockScore").getTextContent()).doubleValue();
		blockScore[blockNr] = thisBlockScore;

		double thisBlockRmsd = new Double(map.getNamedItem("blockRmsd").getTextContent()).doubleValue();
		blockRmsd[blockNr] = thisBlockRmsd;


		// parse out the equivalent positions from the file
		int nrEqr = 0;
		for ( int e =0; e< numChildren ; e++){
			Node  eqr = valList.item(e);

			if(!eqr.hasAttributes()) continue;


			if ( eqr.getNodeName().equals("eqr")) {
				nrEqr++;
				NamedNodeMap atts = eqr.getAttributes();

				int eqrNr = new Integer(atts.getNamedItem("eqrNr").getTextContent()).intValue();

				String pdbres1 = atts.getNamedItem("pdbres1").getTextContent();
				String chain1 = atts.getNamedItem("chain1").getTextContent();
				String pdbres2 = atts.getNamedItem("pdbres2").getTextContent();
				String chain2 = atts.getNamedItem("chain2").getTextContent();

				//System.out.println(blockNr + " " + eqrNr + " " + chain1+" " + pdbres1 + ":" + chain2 + " " + pdbres2);

				pdbAln[blockNr][0][eqrNr] = chain1+":"+pdbres1;
				pdbAln[blockNr][1][eqrNr] = chain2+":"+pdbres2;

				//  A WORK AROUND FOR THE PROBLEM THAT WE DON:T HAVE PDBs LOADED AT THIS TIME...

				/* int pos1 = getPositionForPDBresunm(pdbres1,chain1,ca1);
				int pos2 = getPositionForPDBresunm(pdbres2,chain2,ca2);
				//System.out.println("settion optAln " + blockNr + " " + eqrNr + " " + pos1);
				optAln[blockNr][0][eqrNr] = pos1;
				optAln[blockNr][1][eqrNr] = pos2;
				 */
			} else if ( eqr.getNodeName().equals("matrix")){
				// process Matrix
				Matrix m = new Matrix(3,3);

				for (int i =1 ; i <= 3 ; i++){
					for (int j =1 ; j <= 3 ; j++){
						String att = getAttribute(eqr, "mat" +i + j);
						double val = Double.parseDouble(att);
						m.set(i-1,j-1,val);

					}	
				}
				ms[blockNr] = m;

			} else if ( eqr.getNodeName().equals("shift")){
				Atom shift = new AtomImpl();
				double x = Double.parseDouble(getAttribute(eqr, "x"));
				double y = Double.parseDouble(getAttribute(eqr, "y"));
				double z = Double.parseDouble(getAttribute(eqr, "z"));
				shift.setX(x);
				shift.setY(y);
				shift.setZ(z);
				shifts[blockNr] = shift;

			}

		}
		//System.out.println("setting block " + blockNr + " eqr: " + nrEqr);
		optLen[blockNr] = nrEqr;



		a.setOptLen(optLen);
		//a.setOptAln(optAln);
		a.setPdbAln(pdbAln);
		a.setBlockGap(blockGap);
		a.setBlockSize(blockSize);

		a.setBlockScore(blockScore);
		a.setBlockRmsd(blockRmsd);



	}

	private static String getAttribute(Node node, String attr){
		if( ! node.hasAttributes()) 
			return null;

		NamedNodeMap atts = node.getAttributes();

		if ( atts == null)
			return null;

		Node att = atts.getNamedItem(attr);
		if ( att == null)
			return null;

		String value = att.getTextContent();

		return value;

	}

	/** get the position of PDB residue nr X in the ato marray
	 * 
	 * @param pdbresnum
	 * @param chainId
	 * @param atoms
	 * @return
	 */
	private static int getPositionForPDBresunm(String pdbresnum, String chainId , Atom[] atoms){
		ResidueNumber residueNumber =  ResidueNumber.fromString(pdbresnum);
		residueNumber.setChainId(chainId);
		
		for ( int i =0; i< atoms.length ;i++){
			Group g = atoms[i].getGroup();
			
			//System.out.println(g.getResidueNumber() + "< ? >" + residueNumber +"<");
			if ( g.getResidueNumber().equals(residueNumber)){
				//System.out.println(g + " == " + residueNumber );
				Chain c = g.getChain();
				if ( c.getChainID().equals(chainId)){
					return i;
				}
			}
		}
		return -1;
	}



}

