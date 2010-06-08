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
      

      if ( ! (parameters instanceof SmithWaterman3DParameters))
         throw new IllegalArgumentException("provided parameter object is not of type SmithWaterman3DParameters");
      
      params = (SmithWaterman3DParameters) parameters;
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
               params.getMatch(),
               params.getReplace(),
               params.getInsert(),
               params.getDelete(),
               params.getGapExtend(),
               
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

   /**
    * Converts a sequence alignment into a structural alignment
    * @param aligner The sequence aligner
    * @param ca1 CA atoms from the query sequence
    * @param ca2 CA atoms from the target sequence
    * @param aligPair The sequence alignment calculated by aligner
    * @return an AFPChain encapsulating the alignment in aligPair
    * @throws StructureException
    */
   private AFPChain convert(AlignmentAlgorithm aligner,Atom[] ca1, Atom[] ca2,  AlignmentPair aligPair) throws StructureException {
      AFPChain afpChain = new AFPChain();
      int ca1Length = ca1.length;
      int ca2Length = ca2.length;		

      afpChain.setAlignScore(aligPair.getScore());
      afpChain.setCa1Length(ca1Length);
      afpChain.setCa2Length(ca2Length);

      int nAtom=0; 
      int nGaps=0;

      List<String> labels = aligPair.getLabels();

      if ( ! (labels.size() == 2)){
         throw new StructureException("Expected labels of length 2 but got " + labels.size());
      }

      SymbolList symb1 = aligPair.symbolListForLabel(labels.get(0));
      SymbolList symb2 = aligPair.symbolListForLabel(labels.get(1));
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


      afpChain.setGapLen(nGaps);
      afpChain.setAlnseq1(alnseq1);
      //System.out.println(new String(alnseq1));
      //System.out.println(new String(alnsymb));
      //System.out.println(new String(alnseq2));
      afpChain.setAlnseq2(alnseq2);
      afpChain.setAlnsymb(alnsymb);


      // CE uses the aligned pairs as reference not the whole alignment including gaps...
      afpChain.setIdentity(nrIdent*1.0/pos);
      afpChain.setSimilarity(nrSim*1.0/pos);

      afpChain.setAlnLength(aligLength);
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
      afpChain.setAlnbeg1(queryStart-1);
      afpChain.setAlnbeg2(subjectStart-1);
      cecalc.setAlign_se1(align_se1);
      cecalc.setAlign_se2(align_se2);
      //System.out.println(align_se1);
      //System.out.println(align_se2);

      //System.out.println("lcmp:" + lcmp + " " + aligLength + " " + nAtom);
      cecalc.setLcmp(aligLength);
      
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
