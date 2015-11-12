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
package org.biojava.nbio.structure.align.seq;


import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.AbstractStructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.ce.UserArgumentProcessor;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.ConfigurationException;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** provides a 3D superimposition based on the sequence alignment
 * 
 * @author Andreas Prlic
 *
 */
public class SmithWaterman3Daligner extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "Smith-Waterman superposition";
	
	private static final Logger logger = LoggerFactory.getLogger(SmithWaterman3Daligner.class);

	/**
	 *  version history:
	 *  1.1 - Added more parameters to the command line, including -maxOptRMSD
	 *  1.0 - Initial version
	 */
	private static final String version = "1.1";
	SmithWaterman3DParameters params;

	public static void main(String[] args) throws ConfigurationException {
		//args = new String[]{"-pdb1","1cdg.A","-pdb2","1tim.A","-pdbFilePath","/tmp/","-show3d","-printFatCat"};
		UserArgumentProcessor processor = new SmithWatermanUserArgumentProcessor();
		processor.process(args);
	}

	public SmithWaterman3Daligner(){
		params = new SmithWaterman3DParameters();
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
		if ( params == null)
			params = new SmithWaterman3DParameters();
		return align(ca1,ca2,params);
	}

	@Override
	public AFPChain align(Atom[] ca1, Atom[] ca2, Object parameters)
	throws StructureException {

		if ( parameters == null) {
			throw new IllegalArgumentException("Got null instead of SmithWaterman3DParameters!");
		}
		if ( ! (parameters instanceof SmithWaterman3DParameters))
			throw new IllegalArgumentException("provided parameter object is not of type SmithWaterman3DParameters, but " + parameters.getClass().getName());

		params = (SmithWaterman3DParameters) parameters;
		AFPChain afpChain = new AFPChain();

		
			// covert input to sequences
			String seq1 = StructureTools.convertAtomsToSeq(ca1);
			String seq2 = StructureTools.convertAtomsToSeq(ca2);

			ProteinSequence s1 = null;
			ProteinSequence s2 = null;
					
			try {
				s1 = new ProteinSequence(seq1);
				s2 = new ProteinSequence(seq2);
			} catch (CompoundNotFoundException e){

				throw new StructureException(e.getMessage(),e);
			}

			// default blosum62 
			SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();
					
			GapPenalty penalty = new SimpleGapPenalty();
			
			penalty.setOpenPenalty(params.getGapOpen());
			penalty.setExtensionPenalty(params.getGapExtend());
			
			PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

			SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();
			
			if (pair.getTarget().toString().isEmpty() || pair.getQuery().toString().isEmpty()) {
				throw new StructureException("Empty alignment for sequences "+s1+" and "+s2);
			}

			logger.debug("Smith-Waterman alignment is: "+pair.toString(100));
			
			// convert to a 3D alignment...
			afpChain = convert(ca1,ca2,pair, smithWaterman);
			


		return afpChain;
	}

	/**
	 * Converts a sequence alignment into a structural alignment
	 * @param smithWaterman The sequence aligner
	 * @param ca1 CA atoms from the query sequence
	 * @param ca2 CA atoms from the target sequence
	 * @param smithWaterman pairwise Sequence aligner
	 * @param pair The sequence alignment calculated by aligner
	 * @return an AFPChain encapsulating the alignment in aligPair
	 * @throws StructureException
	 */
	private AFPChain convert(Atom[] ca1, Atom[] ca2,  SequencePair<ProteinSequence, 
			AminoAcidCompound> pair, PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman) throws StructureException {
		AFPChain afpChain = new AFPChain();
		int ca1Length = ca1.length;
		int ca2Length = ca2.length;		

		afpChain.setAlignScore(smithWaterman.getScore());
		afpChain.setCa1Length(ca1Length);
		afpChain.setCa2Length(ca2Length);

		int nrCols = pair.getLength();
		int nAtom=0; 
		int nGaps=0;

		Atom[] strBuf1 = new Atom[nrCols];
		Atom[] strBuf2 = new Atom[nrCols];
		char[] alnseq1 = new char[ca1Length+ca2Length+1];
		char[] alnseq2 = new char[ca1Length+ca2Length+1] ;
		char[] alnsymb = new char[ca1Length+ca2Length+1];

		Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");

		int pos = 0 ; // aligned positions
		
		int nrIdent = 0;
		int nrSim   = 0;

		int[] align_se1 = new int[nrCols+1];
		int[] align_se2 = new int[nrCols+1];

		for ( int i = 1 ; i <= nrCols; i ++){

			int myI = i-1;

			Compound s1 =  pair.getCompoundAt(1, i);
			Compound s2 =  pair.getCompoundAt(2,i);

			// alignment is using internal index start at 1...
			int pos1 = pair.getIndexInQueryAt(i)  -1;
			int pos2 = pair.getIndexInTargetAt(i) -1;
			
			if ( ( ! s1.equals(gapSymbol) )&&  (! s2.equals(gapSymbol))){
				strBuf1[nAtom] = ca1[pos1];
				strBuf2[nAtom] = ca2[pos2];
				//
				char l1 = getOneLetter(ca1[pos1].getGroup());
				char l2 = getOneLetter(ca2[pos2].getGroup());
				//
				alnseq1[myI] = Character.toUpperCase(l1);
				alnseq2[myI] = Character.toUpperCase(l2);
				alnsymb[myI] = ' ';
				//
				if ( l1 == l2) {					
					nrIdent++;
					nrSim++;
					alnsymb[myI] = '|';

				} else if ( AFPAlignmentDisplay.aaScore(l1, l2) > 0){

					nrSim++;
					alnsymb[myI] = ':';
				}
				//
				align_se1[myI] = pos1;
				align_se2[myI] = pos2;
				//
				pos++;
				nAtom++;

			} else {
				// there is a gap at this position
				nGaps++;

				alnsymb[myI] = ' ';
				align_se1[myI] = -1;
				align_se2[myI] = -1;

				if ( s1.equals(gapSymbol)){
					alnseq1[myI] = '-';

				} else {
					char l1 = getOneLetter(ca1[pos1].getGroup());
					alnseq1[myI] = Character.toUpperCase(l1);
					align_se1[myI] = pos1;
				}
				if ( s2.equals(gapSymbol)){
					alnseq2[myI] = '-';

				} else {
					char l2 = getOneLetter(ca2[pos2].getGroup());
					alnseq2[myI] = Character.toUpperCase(l2);
					align_se2[myI] = pos2;

				}
			}

		}

		afpChain.setAlnbeg1(pair.getIndexInQueryAt(1)-1);
		afpChain.setAlnbeg2(pair.getIndexInTargetAt(1)-1);
		
		afpChain.setGapLen(nGaps);
		afpChain.setAlnseq1(alnseq1);

		afpChain.setAlnseq2(alnseq2);
		afpChain.setAlnsymb(alnsymb);

		// CE uses the aligned pairs as reference not the whole alignment including gaps...
		afpChain.setIdentity(nrIdent*1.0/pos);
		afpChain.setSimilarity(nrSim*1.0/pos);

		afpChain.setAlnLength(nrCols);
		afpChain.setOptLength(nAtom);
		int[] optLen = new int[]{nAtom};
		afpChain.setOptLen(optLen);
		
		if(nAtom<4) 
			return afpChain;

		CeParameters params = new CeParameters();
		CECalculator cecalc = new CECalculator(params);

		// here we don't store the rotation matrix for the user!

		double rmsd= cecalc.calc_rmsd(strBuf1, strBuf2, nAtom,true);

		afpChain.setBlockRmsd(new double[]{rmsd});
		afpChain.setOptRmsd(new double[]{rmsd});
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setChainRmsd(rmsd);

		// let's hijack the CE implementation
		// and use some utilities from there to 
		// build up the afpChain object

		cecalc.setnAtom(nAtom);

		cecalc.setAlign_se1(align_se1);
		cecalc.setAlign_se2(align_se2);

		cecalc.setLcmp(nrCols  );

		cecalc.convertAfpChain(afpChain, ca1, ca2);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		return afpChain;
	}

	private static char getOneLetter(Group g){

		if (g==null) return StructureTools.UNKNOWN_GROUP_LABEL;
		
		return StructureTools.get1LetterCode(g.getPDBName());
		
	}

	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		return params;
	}

	@Override
	public String getVersion() {
		return version;
	}

	@Override
	public void setParameters(ConfigStrucAligParams parameters) {
		if ( ! (parameters instanceof SmithWaterman3DParameters))
			throw new IllegalArgumentException("provided parameter object is not of type SmithWaterman3DParameters");
		params = (SmithWaterman3DParameters)parameters;
	}



}
