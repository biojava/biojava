import org.biojava3.aaproperties.IPeptideProperties;
import org.biojava3.aaproperties.PeptidePropertiesImpl;
import org.biojava3.core.sequence.ProteinSequence;

/**
 * Created by andreas on 8/9/14.
 */


public class BioJavaAADemo {

    public static void main(String[] args){
        ProteinSequence pSequence = new ProteinSequence("VLSPADKTNVKAAWGKVGAHAG");

        IPeptideProperties pp = new PeptidePropertiesImpl();

        System.out.println(pp.getIsoelectricPoint(pSequence));;
    }
}
