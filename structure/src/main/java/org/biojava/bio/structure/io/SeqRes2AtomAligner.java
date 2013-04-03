/*
 *                  BioJava development code
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
 * Created on Aug 23, 2007
 *
 */

package org.biojava.bio.structure.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;



/** Aligns the SEQRES residues to the ATOM residues.
 * The AminoAcids that can be matched between the two of them will be set in the SEQRES
 * chains
 *
 *
 * @author Andreas Prlic
 *
 */
public class SeqRes2AtomAligner {

    boolean DEBUG = false;

    static final List<String> excludeTypes;
    private static final FiniteAlphabet alphabet;
    private static final Symbol gapSymbol ;
    private static  SubstitutionMatrix matrix;


    String alignmentString;
    static {
        excludeTypes = new ArrayList<String>();
        excludeTypes.add("HOH"); // we don't want to align water against the SEQRES...
        excludeTypes.add("DOD"); // deuterated water


        alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
        gapSymbol =  alphabet.getGapSymbol();

        try {
            matrix = getSubstitutionMatrix(alphabet);

        } catch (BioException e){
            e.printStackTrace();
        } catch (IOException e){
            e.printStackTrace();
        }
        //matrix.printMatrix();

    }



    public SeqRes2AtomAligner(){
        alignmentString = "";
    }

    public String getAlignmentString() {
        return alignmentString;
    }

    public boolean isDEBUG() {
        return DEBUG;
    }

    public void setDEBUG(boolean debug) {
        DEBUG = debug;
    }

    private Chain getMatchingAtomRes(Chain seqRes, List<Chain> atomList)
    throws StructureException {
        Iterator<Chain> iter = atomList.iterator();
        while(iter.hasNext()){
            Chain atom = iter.next();
            if ( atom.getName().equals(seqRes.getName())){
                return atom;
            }
        }
        throw new StructureException("could not match seqres chainID >" + seqRes.getName() + "< to ATOM chains!");
    }
    public void align(Structure s, List<Chain> seqResList){

        //List<Chain> seqResList = s.getSeqRes();
        List<Chain> atomList   = s.getModel(0);

        Iterator<Chain> iter = seqResList.iterator();
       // List<Chain> chains = new ArrayList<Chain>();
        while ( iter.hasNext()){
            Chain seqRes = iter.next();

            if ( seqRes.getAtomGroups("amino").size() < 1) {
            	if (DEBUG){
            		System.out.println("chain " + seqRes.getName() + " does not contain amino acids, ignoring...");
            	}
            	continue;
            }

            try {

                Chain atomRes = getMatchingAtomRes(seqRes,atomList);
                if ( atomRes.getAtomGroups("amino").size() < 1) {
                	if (DEBUG){
                		System.out.println("chain " + atomRes.getName() + " does not contain amino acids, ignoring...");
                	}
                	continue;
                }
                if ( DEBUG )
                    System.out.println("Alignment for chain "+ atomRes.getName());

                List<Group> seqResGroups = seqRes.getAtomGroups();
                boolean noMatchFound = align(seqResGroups,atomRes.getAtomGroups());
                if ( ! noMatchFound){
                    atomRes.setSeqResGroups(seqResGroups);
                }

                //chains.add(mapped);
            } catch (StructureException e){
                e.printStackTrace();
            }
        }
        //s.setChains(0,chains);



    }


    /** returns the full sequence of the Atom records of a chain
     * with X instead of HETATMSs. The advantage of this is
     * that it allows us to also align HETATM groups back to the SEQRES.
     * @param groups the list of groups in a chain
     *
     * @return string representations
     */
    public String getFullAtomSequence(List<Group> groups){


        StringBuffer sequence = new StringBuffer() ;
        for ( int i=0 ; i< groups.size(); i++){
            Group g = (Group) groups.get(i);
            if ( g instanceof AminoAcid ){
                AminoAcid a = (AminoAcid)g;
                sequence.append( a.getAminoType());
            } else {
                if ( ! excludeTypes.contains(g.getPDBName()))
                    sequence.append("X");
            }

        }

        return sequence.toString();

    }

    /** aligns two chains of groups, where the first chain is representing the
     * list of amino acids as obtained from the SEQRES records, and the second chain
     * represents the groups obtained from the ATOM records (and containing the actual ATOM information).
     * This does the actual alignment and if a group can be mapped to a position in the SEQRES then the corresponding
     * position is repplaced with the group that contains the atoms.
     *
     * @param seqRes
     * @param atomRes
     * @return true if no match has bee found
     * @throws StructureException
     */
    public boolean align(List<Group> seqRes, List<Group> atomRes) throws StructureException{
       /** int MAX_SIZE = 1000;
        if ( (seqRes.size() > MAX_SIZE)
                ||( atomRes.size() > MAX_SIZE) ) {
                    System.err.println("can not align chains, length size exceeds limits!");
                    return false;
                }
        */
        String seq1 = getFullAtomSequence(seqRes);
        //String seq1 = seqRes.getSeqResSequence();
        String seq2 = getFullAtomSequence(atomRes);

        //System.out.println("align seq1 " + seq1);
        //System.out.println("align seq2 " + seq2);


        SimpleAlignment simpleAli = null;
        try  {


            Sequence bjseq1 = ProteinTools.createProteinSequence(seq1,"seq1");
            Sequence bjseq2 = ProteinTools.createProteinSequence(seq2,"seq2");

            //System.out.println(bjseq1.getAlphabet());
            NeedlemanWunsch aligner = new NeedlemanWunsch((short)-2,(short) 5,(short) 3,(short) 3,(short) 0,matrix);

            if (DEBUG){
            	System.out.println("seq lengths: " + bjseq1.seqString().length() + " " + bjseq2.seqString().length());
            	System.out.println("seq1: " + bjseq1.seqString());
            	System.out.println("seq2: " + bjseq2.seqString());
            }
         
            org.biojava.bio.alignment.Alignment ali = aligner.getAlignment(bjseq1,bjseq2);
            if ( ! (ali instanceof SimpleAlignment )) {
                throw new Exception ("Alignment is not a SimpleAlignment!");

            }

            simpleAli = (SimpleAlignment) ali;
     
            alignmentString = aligner.getAlignmentString();
     
            if (DEBUG)
                System.out.println(alignmentString);
        } catch (Exception e){
            e.printStackTrace();
            System.err.println("align seq1 " + seq1);
            System.err.println("align seq2 " + seq2);

        }

        if ( simpleAli == null)
            throw new StructureException("could not align objects!");

        //System.out.println(ali.getAlphabet());
        SymbolList lst1 = simpleAli.symbolListForLabel("seq1");
        //System.out.println("from seqres : " + lst1.seqString());
        SymbolList lst2 = simpleAli.symbolListForLabel("seq2");

        boolean noMatchFound = mapChains(seqRes,lst1,atomRes,lst2,gapSymbol);
        return noMatchFound;

    }


    private boolean mapChains(List<Group> seqRes, SymbolList lst1,
            List<Group> atomRes, SymbolList lst2,Symbol gapSymbol) throws StructureException{

        assert (lst1.length() == lst2.length());

        // at the present stage the future seqREs are still stored as Atom groups in the seqRes chain...
        List<Group> seqResGroups = seqRes;

        int aligLength = lst1.length();
        int posSeq  = -1;
        int posAtom = -1;
        // System.out.println(gapSymbol.getName());

        // make sure we actually find an alignment
        boolean noMatchFound = true;

        for (int i = 1; i <= aligLength; i++) {

            Symbol s = lst1.symbolAt(i);
            Symbol a = lst2.symbolAt(i);

            //TODO replace the text gap with the proper symbol for terminal gap

            // don't count gaps and terminal gaps
            String sn = s.getName();
            String an = a.getName();
            if (! ( sn.equals("gap")  || sn.equals(gapSymbol.getName()))){
                posSeq++;
            }
            if (! ( an.equals("gap") || an.equals(gapSymbol.getName()))){
                posAtom++;
            }


           /* System.out.println(i+ " " + posSeq + " " +
                    s.getName() +
                    " " + posAtom +
                    " " + a.getName() +
                    " " + s.getName().equals(a.getName()));
           */

            if ( s.getName().equals(a.getName())){
                // the atom record can be aligned to the SeqRes record!
                // replace the SeqRes group with the Atom group!

                Group s1 = seqRes.get(posSeq);
                Group a1 = atomRes.get(posAtom);
                //System.out.println(s1.getPDBName() + " == " + a1.getPDBName());
                // need to trim the names to allow matching e.g in
                // pdb1b2m
                String pdbNameS = s1.getPDBName();
                String pdbNameA = a1.getPDBName();
                if ( pdbNameS == null || pdbNameA == null ){
                	System.err.println("nullvalue found at " + posSeq + " when trying to align " + s1 + " and " + a1 + " " + posAtom);
                	throw new StructureException("nullvalue found at group.getPDBName()");
                }
                if (! pdbNameA.equals(pdbNameS)){
                    if ( ! pdbNameA.trim().equals(pdbNameS.trim())) {
                        System.err.println(s1 + " " + posSeq + " does not align with " + a1+ " " + posAtom);
                        //System.exit(0);// for debug only
                        throw new StructureException("could not match residues");

                    }
                }

                // do the actual replacing of the SEQRES group with the ATOM group
                seqResGroups.set(posSeq,a1);
                noMatchFound = false;
            }
        }


        // now we merge the two chains into one
        // the Groups that can be aligned are now pointing to the
        // groups in the Atom records.
        if (  noMatchFound) {

            if ( DEBUG )
                System.out.println("no alignment found!");
        }
        return noMatchFound;

    }

    public static SubstitutionMatrix getSubstitutionMatrix(FiniteAlphabet alphabet)
    throws IOException,BioException {

        InputStream inStream = SeqRes2AtomAligner.class.getResourceAsStream("/org/biojava/bio/structure/blosum62.mat");

        String newline = System.getProperty("line.separator");
        BufferedReader reader = new BufferedReader(new InputStreamReader( inStream));
        StringBuffer file = new StringBuffer();
        while (reader.ready()){
            file.append(reader.readLine() );
            file.append(newline);
        }


        SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet,file.toString(),"blosum 62");
        return matrix;


    }

}
