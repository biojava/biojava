import org.biojava3.aaproperties.IPeptideProperties;
import org.biojava3.aaproperties.PeptidePropertiesImpl;
import org.biojava3.core.sequence.ProteinSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Created by andreas on 8/9/14.
 */


public class BioJavaAADemo {

	private final static Logger logger = LoggerFactory.getLogger(BioJavaAADemo.class);

    public static void main(String[] args){
        ProteinSequence pSequence = new ProteinSequence("VLSPADKTNVKAAWGKVGAHAG");

        IPeptideProperties pp = new PeptidePropertiesImpl();

        logger.info("Peptide Properties: {}", pp.getIsoelectricPoint(pSequence));;
    }
}
