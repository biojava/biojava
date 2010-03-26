package org.biojava.bio.structure.align.seq;


import java.io.IOException;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioException;

import org.biojava.bio.alignment.AlignmentAlgorithm;
import org.biojava.bio.alignment.AlignmentPair;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

import org.biojava.bio.structure.align.ce.UserArgumentProcessor;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.io.SeqRes2AtomAligner;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;


/** provides a 3D superimposition based on the sequence alignment
 * 
 * @author Andreas Prlic
 *
 */
public class SmithWaterman3Daligner extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "Smith-Waterman superposition";

	private static final String version = "1.0";


	public static void main(String[] args){
		SmithWaterman3Daligner algorithm = new SmithWaterman3Daligner();
		if (args.length  == 0 ) {			
			System.out.println(algorithm.printHelp());
			return;			
		}

		if ( args.length == 1){
			if (args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("-help")|| args[0].equalsIgnoreCase("--help")){
				System.out.println(algorithm.printHelp());								
				return;
			}			
		}

		UserArgumentProcessor processor = new SmithWatermanUserArgumentProcessor();
		processor.process(args);
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
		return align(ca1,ca2,null);
	}



	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object params)
	throws StructureException {
		AFPChain afpChain = new AFPChain();
		try {
			// covert input to sequences
			String seq1 = StructureTools.convertAtomsToSeq(ca1);
			String seq2 = StructureTools.convertAtomsToSeq(ca2);

			// The alphabet of the sequences. 
			//FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");

			// read blosum62 
			SubstitutionMatrix matrix = getBlosum62();


			// Perform a local alginment from the sequences with Smith-Waterman. 
			// Firstly, define the expenses (penalties) for every single operation.
			SmithWaterman aligner = new SmithWaterman(
					(short) -1,     // match
					(short) 3,      // replace 
					(short) 2,      // insert
					(short) 2,      // delete
					(short) 1,      // gapExtend
					matrix  // SubstitutionMatrix
			);

			Sequence query  = ProteinTools.createProteinSequence(seq1,  "seq1");
			Sequence target = ProteinTools.createProteinSequence(seq2,  "seq2");
			
			AlignmentPair aligPair = aligner.pairwiseAlignment(query, target);

			// Print the alignment to the screen
			// System.out.println(aligner.getAlignmentString());
			// convert to a 3D alignment...
			 afpChain = convert(aligner,ca1,ca2,aligPair);
			
		} catch (Exception e){

			throw new StructureException(e.getMessage(),e);
		}
		return afpChain;
	}

	private AFPChain convert(AlignmentAlgorithm aligner,Atom[] ca1, Atom[] ca2,  AlignmentPair alig) throws StructureException {
		AFPChain afpChain = new AFPChain();
		int ca1Length = ca1.length;
		int ca2Length = ca2.length;		
				
		afpChain.setAlignScore(alig.getScore());

		alig.getLabels();

		int nAtom=0; 
		int nGaps=0;
		
		List<String> labels = alig.getLabels();
		
		if ( ! (labels.size() == 2)){
			throw new StructureException("Expected labels of length 2 but got " + labels.size());
		}
		SymbolList symb1 = alig.symbolListForLabel(labels.get(0));
		SymbolList symb2 = alig.symbolListForLabel(labels.get(1));
		
		List<Alphabet> alphas = new ArrayList<Alphabet>();
		alphas.add(symb1.getAlphabet()); 
		Symbol gapSymbol = AlphabetManager.getGapSymbol(alphas);
		
		int lcmp = symb1.length();
		
		Atom[] strBuf1 = new Atom[lcmp];
		Atom[] strBuf2 = new Atom[lcmp];
		
		int pos1 = -1;
		int pos2 = -1;
		
		char[] alnseq1 = new char[ca1Length+ca2Length+1];
		char[] alnseq2 = new char[ca1Length+ca2Length+1] ;
		char[] alnsymb = new char[ca1Length+ca2Length+1];
		
		int nrIdent = 0;
		int nrSim   = 0;
		int pos =0;
		
		int[] align_se1 = new int[lcmp];
		int[] align_se2 = new int[lcmp];
		
		for(int ia=0; ia<lcmp; ia++) {
			Symbol s1 = symb1.symbolAt(ia+1);
			Symbol s2 = symb2.symbolAt(ia+1);
			if ( !s1.equals(gapSymbol))
				pos1++;
			if ( ! s2.equals(gapSymbol))
				pos2++;
			
			if ( ( ! s1.equals(gapSymbol) )&&  (! s2.equals(gapSymbol))){
				
				strBuf1[nAtom]=ca1[pos1];
				strBuf2[nAtom]=ca2[pos2];
				
				char l1 = getOneLetter(ca1[pos1].getParent());
				char l2 = getOneLetter(ca2[pos2].getParent());

				alnseq1[ia] = Character.toUpperCase(l1);
				alnseq2[ia] = Character.toUpperCase(l2);
				alnsymb[ia] = ' ';
				if ( l1 == l2) {					
					nrIdent++;
					nrSim++;
					alnsymb[ia] = '|';
				} else if ( AFPAlignmentDisplay.aaScore(l1, l2) > 0){
					nrSim++;
					alnsymb[ia] = ':';
				}
		
				align_se1[ia] = pos1;
				align_se2[ia] = pos2;
				pos++;
				nAtom++;
			} else {
				// there is a gap at this position
				nGaps++;

				alnsymb[ia] = ' ';
				align_se1[ia] = -1;
				align_se2[ia] = -1;
				
				if ( s1.equals(gapSymbol)){
					alnseq1[ia] = '-';
					
				} else {
					char l1 = getOneLetter(ca1[pos1].getParent());
					alnseq1[ia] = Character.toUpperCase(l1);
					align_se1[ia] = pos1;
				}
				if ( s2.equals(gapSymbol)){
					alnseq2[ia] = '-';
					
				} else {
					char l2 = getOneLetter(ca2[pos2].getParent());
					alnseq2[ia] = Character.toUpperCase(l2);
					align_se2[ia] = pos2;
					
				}
				
			}
			
		}

		
		afpChain.setGapLen(nGaps);
		afpChain.setAlnseq1(alnseq1);
		afpChain.setAlnseq2(alnseq2);
		afpChain.setAlnsymb(alnsymb);
		
		
		// CE uses the aligned pairs as reference not the whole alignment including gaps...
		afpChain.setIdentity(nrIdent*1.0/pos);
		afpChain.setSimilarity(nrSim*1.0/pos);
		
		afpChain.setAlnLength(nAtom);
		afpChain.setOptLength(nAtom);
		int[] optLen = new int[]{nAtom};
		afpChain.setOptLen(optLen);
		
		if(nAtom<4) 
			return afpChain;
		
		CeParameters params = new CeParameters();
		CECalculator cecalc = new CECalculator(params);
		//sup_str(strBuf1, strBuf2, nAtom, _d);
		// here we don't store the rotation matrix for the user!
		double rmsd= cecalc.calc_rmsd(strBuf1, strBuf2, nAtom,true, false);
		
		afpChain.setBlockRmsd(new double[]{rmsd});
		afpChain.setOptRmsd(new double[]{rmsd});
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setChainRmsd(rmsd);
	
		
		// let's hijack the CE implementation
		// and use some utilities from there to 
		// build up the afpChain object
		
		cecalc.setnAtom(nAtom);
		cecalc.setLcmp(lcmp);
		cecalc.setAlign_se1(align_se1);
		cecalc.setAlign_se2(align_se2);
		
		cecalc.convertAfpChain(afpChain, ca1, ca2);
				
		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);
		//AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2);
		
		return afpChain;
	}
	
	private static char getOneLetter(Group g){

		try {
			Character c = StructureTools.get1LetterCode(g.getPDBName());
			return c;
		} catch (Exception e){
			return 'X';
		}
	}

	private SubstitutionMatrix getBlosum62() throws BioException,IOException {

		return SeqRes2AtomAligner.getSubstitutionMatrix((FiniteAlphabet)AlphabetManager.alphabetForName("PROTEIN-TERM"));

	}

	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getVersion() {
		return version;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		// no parameters as of yet.
		
	}



}
