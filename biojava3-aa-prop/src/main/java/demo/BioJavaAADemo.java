package demo;
import org.biojava3.aaproperties.IPeptideProperties;
import org.biojava3.aaproperties.PeptidePropertiesImpl;
import org.biojava3.core.exceptions.CompoundNotFoundException;
import org.biojava3.core.sequence.ProteinSequence;

/**
 * Created by andreas on 8/9/14.
 */


public class BioJavaAADemo {


    public static void main(String[] args) throws CompoundNotFoundException {
        ProteinSequence pSequence = new ProteinSequence("VLSPADKTNVKAAWGKVGAHAG");

        IPeptideProperties pp = new PeptidePropertiesImpl();

        System.out.println("Peptide Properties: " + pp.getIsoelectricPoint(pSequence));
    }
}
