/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.bio.alignment;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 *
 * @author Chris Friedline <cfriedline@vcu.edu>
 */
public class BioJavaAlignmentDemo {

    private final short GAPOPEN = 10;
    private final short GAPEXTEND = 5;
    private final short MATCH = 0;
    private final short MISMATCH = 0;

    public static void main(String[] args) {
        BioJavaAlignmentDemo test = new BioJavaAlignmentDemo();
        test.testrun();
    }

    public void testrun() {
        try {
            String header1 = ">gb|CP001821.1|:690821-692026 translation elongation factor Tu [Xylanimonas cellulosilytica DSM 15894]";
            String protein1 = "VAKAKFERTKPHVNVGTIGHVDHGKTTLTAAISKTLAEKYPASEGYLANQVVDFDGIDKAPEEKQRGITINISHIEYETPNRHYAHVDAPGHADYIKNMITGAAQMDGAILVVAATDGPMAQTREHVLLARQVGVPYLLVALNKSDMVDDEEILELVEMEVRELLSSQGFDGDDAPVVRVSGLKALEGDPEWQAKVLELMEAVDTNVPEPVRDLDKPFLMPIEDVFTITGRGTVVTGKVERGALNVNSEVEIVGIRNPQKTTVTGIETFHKSMDQAQAGDNTGLLLRGIKREDVERGQVVVKPGSITPHTDFEAQVYILGKDEGGRHNPFYSNYRPQFYFRTTDVTGVISLPEGTEMVMPGDNTEMTVELIQPIAMEEGLGFAIREGGRTVGSGRVTKIIK";
            String header2 = ">dbj|AP006618.1|:c5375075-5373828 putative translation elongation factor TU [Nocardia farcinica IFM 10152]";
            String protein2 = "MTPRTAATAGTNTVQEDKTVAKAKFERTKPHVNIGTIGHVDHGKTTLTAAITKVLADKYPDLNQSFAFDQIDKAPEEKARGITINISHVEYQTEKRHYAHVDAPGHADYIKNMITGAAQMDGAILVVAATDGPMPQTREHVLLARQVGVPYILVALNKADMVDDEEILELVEMEVRELLAAQEFDEEAPVVRVSGLKALEGDPKWVKSVEDLMDAVDESIPDPVRETDKPFLMPIEDVFTITGRGTVVTGRVERGIINVNEEVEITGIRPETTKTTVTGIEMFRKLLDQGQAGDNVGLLIRGIKREDVERGQVVIKPGTTTPHTEFEGQAYILSKDEGGRHTPFFNNYRPQFYFRTTDVTGVVTLPEGTEMVMPGDNTEMSVKLIQPVAMEEGLRFAIREGGRTVGAGRVTKIIK";
            Sequence query = ProteinTools.createProteinSequence(protein1, header1.split("\\|")[1]);
            Sequence target = ProteinTools.createProteinSequence(protein2, header2.split("\\|")[1]);
            FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
            SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(getClass().getResource("BLOSUM50").toURI()));
            NeedlemanWunsch aligner = new NeedlemanWunsch(MATCH, MISMATCH, GAPOPEN, GAPOPEN, GAPEXTEND, matrix);
//            Alignment alignment = aligner.getAlignment(query, target);
//            System.out.println(aligner.getAlignmentString());
            AlignmentPair pair = aligner.pairwiseAlignment(query, target);
			System.out.println(pair.formatOutput(100));
            System.out.printf("\n%d\t%d\n", protein1.length(), protein2.length());
        } catch (Exception ex) {
            Logger.getLogger(BioJavaAlignmentDemo.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}
