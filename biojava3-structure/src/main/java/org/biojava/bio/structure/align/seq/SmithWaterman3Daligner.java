package org.biojava.bio.structure.align.seq;


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
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.Compound;


/** provides a 3D superimposition based on the sequence alignment
 * 
 * @author Andreas Prlic
 *
 */
public class SmithWaterman3Daligner extends AbstractStructureAlignment implements StructureAlignment {

	public static final String algorithmName = "Smith-Waterman superposition";

	private static final String version = "1.0";
	SmithWaterman3DParameters params;

	public static void main(String[] args){

		args = new String[]{"-pdb1","1cdg.A","-pdb2","1tim.A","-pdbFilePath","/tmp/","-show3d","-printFatCat"};
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

		try {
			// covert input to sequences
			String seq1 = StructureTools.convertAtomsToSeq(ca1);
			String seq2 = StructureTools.convertAtomsToSeq(ca2);

			//System.out.println(seq1);
			//System.out.println(seq2);
			
			ProteinSequence s1 = new ProteinSequence(seq1);
			ProteinSequence s2 = new ProteinSequence(seq2);

			// default blosum62 
			SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();
					
			GapPenalty penalty = new SimpleGapPenalty();
			
			
			penalty.setOpenPenalty(params.getGapOpen());
			penalty.setExtensionPenalty(params.getGapExtend());
			
			PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
				Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

			SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();

			// Print the alignment to the screen
			// System.out.println(aligner.getAlignmentString());
			// convert to a 3D alignment...
			afpChain = convert(ca1,ca2,pair, smithWaterman);

		} catch (Exception e){

			throw new StructureException(e.getMessage(),e);
		}
		return afpChain;
	}

	/**
	 * Converts a sequence alignment into a structural alignment
	 * @param smithWaterman The sequence aligner
	 * @param ca1 CA atoms from the query sequence
	 * @param ca2 CA atoms from the target sequence
	 * @param smithWaterman pairwise Sequence aligner
	 * @param aligPair The sequence alignment calculated by aligner
	 * @return an AFPChain encapsulating the alignment in aligPair
	 * @throws StructureException
	 */
	private AFPChain convert(Atom[] ca1, Atom[] ca2,  SequencePair<ProteinSequence, 
			AminoAcidCompound> pair, PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman) throws StructureException {
		AFPChain afpChain = new AFPChain();
		int ca1Length = ca1.length;
		int ca2Length = ca2.length;		

		//System.out.println(aligner.getScore());

		afpChain.setAlignScore(smithWaterman.getScore());
		afpChain.setCa1Length(ca1Length);
		afpChain.setCa2Length(ca2Length);

		//int nrRows = pair.getSize();
		int nrCols = pair.getLength();
		int nAtom=0; 
		int nGaps=0;

		//System.out.println(pair.toString(60));

		Atom[] strBuf1 = new Atom[nrCols];
		Atom[] strBuf2 = new Atom[nrCols];
		char[] alnseq1 = new char[ca1Length+ca2Length+1];
		char[] alnseq2 = new char[ca1Length+ca2Length+1] ;
		char[] alnsymb = new char[ca1Length+ca2Length+1];

//		System.out.println("nrRows: " + nrRows);
//		System.out.println("nrCols: " + nrCols) ;

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
				//System.out.println(i+ " " + 1 + " " + 1 + " " + s1 + " " + s2 + " " + pos1 + " " + ca1.length);
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

//		System.out.println("nr aligned positions:" + pos + " " + nAtom);
//		System.out.println(new String(alnseq1));
//		System.out.println(new String(alnsymb));
//		System.out.println(new String(alnseq2));

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
		//sup_str(strBuf1, strBuf2, nAtom, _d);
		// here we don't store the rotation matrix for the user!
		//System.out.println(strBuf1.length + " " + aligLength);
		double rmsd= cecalc.calc_rmsd(strBuf1, strBuf2, nAtom,true, false);

		afpChain.setBlockRmsd(new double[]{rmsd});
		afpChain.setOptRmsd(new double[]{rmsd});
		afpChain.setTotalRmsdOpt(rmsd);
		afpChain.setChainRmsd(rmsd);
		

		

		// let's hijack the CE implementation
		// and use some utilities from there to 
		// build up the afpChain object

		cecalc.setnAtom(nAtom);
		//afpChain.setAlnbeg1(0);
		//afpChain.setAlnbeg2(0);
		cecalc.setAlign_se1(align_se1);
		cecalc.setAlign_se2(align_se2);
		//System.out.println(align_se1);
		//System.out.println(align_se2);

		//System.out.println("lcmp:" + lcmp + " " + aligLength + " " + nAtom);
		cecalc.setLcmp(nrCols  );

		cecalc.convertAfpChain(afpChain, ca1, ca2);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		return afpChain;
	}



	/*
      //Sequence symb1 = aligPair.getQuery();
      //Sequence symb2 = aligPair.getSubject();

      int queryEnd = aligPair.getQueryEnd();
      int queryStart = aligPair.getQueryStart();
      int subjectEnd = aligPair.getSubjectEnd();
      int subjectStart = aligPair.getSubjectStart();

//      try {
//         System.out.println(aligPair.formatOutput(89));
//      } catch (Exception e) {
//         e.printStackTrace();
//         // TODO: handle exception
//      }
      List<Alphabet> alphas = new ArrayList<Alphabet>();
      alphas.add(symb1.getAlphabet()); 
      Symbol gapSymbol = AlphabetManager.getGapSymbol(alphas);

      //System.out.println(aligPair.getQueryStart() + " " + " " + aligPair.getQueryEnd()+ " "+ aligPair.getQueryLength() + " " + aligPair.getSubjectStart() + " " + aligPair.getSubjectEnd() + " " + aligPair.getSubjectLength() );

      int lcmp = aligPair.getQueryLength() - 1;

      Atom[] strBuf1 = new Atom[lcmp];
      Atom[] strBuf2 = new Atom[lcmp];

      int pos1 = queryStart - 2 ;
      int pos2 = subjectStart - 2;

      char[] alnseq1 = new char[ca1Length+ca2Length+1];
      char[] alnseq2 = new char[ca1Length+ca2Length+1] ;
      char[] alnsymb = new char[ca1Length+ca2Length+1];

      int nrIdent = 0;
      int nrSim   = 0;
      int pos =0;

      int[] align_se1 = new int[lcmp];
      int[] align_se2 = new int[lcmp];
      int aligLength = 0;

      for(int ia=0;ia < Math.min(queryEnd - queryStart, subjectEnd
            - subjectStart) ; ia++) {

         aligLength++;
         Symbol s1 = symb1.symbolAt(ia+queryStart);
         Symbol s2 = symb2.symbolAt(ia+subjectStart);



         if ( !s1.equals(gapSymbol))
            pos1++;
         if ( ! s2.equals(gapSymbol))
            pos2++;

         if ( ( ! s1.equals(gapSymbol) )&&  (! s2.equals(gapSymbol))){
            //System.out.println(ia+ " " + queryStart + " " + subjectStart + " " + s1.getName() + " " + s2.getName() + pos1 + " " + ca1.length);
            strBuf1[nAtom] = ca1[pos1];
            strBuf2[nAtom] = ca2[pos2];

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


      //AFPAlignmentDisplay.getAlign(afpChain, ca1, ca2);
	 */



	private static char getOneLetter(Group g){

		try {
			Character c = StructureTools.get1LetterCode(g.getPDBName());
			return c;
		} catch (Exception e){
			return 'X';
		}
	}



	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public ConfigStrucAligParams getParameters() {
		// TODO Auto-generated method stub
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
