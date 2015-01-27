package org.biojava3.alignment.aaindex;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** The default provider for AAINDEX loads substitution matrices from the AAINDEX file in the resources directory
 * 
 * @author Andreas Prlic
 *
 */
public class DefaultAAIndexProvider implements AAIndexProvider {

	private final static Logger logger = LoggerFactory.getLogger(DefaultAAIndexProvider.class);

	Map<String,SubstitutionMatrix<AminoAcidCompound>> matrices;
	
	public DefaultAAIndexProvider(){
		
		
		InputStream inStream = getInputStreamToAAindexFile();
		
		AAIndexFileParser parser = new AAIndexFileParser();
		
		try {
			parser.parse(inStream);
		} catch (IOException e){
			logger.error("Exception: ", e);
		}
		
		matrices = parser.getMatrices();
		
	}
	
	@Override
	public SubstitutionMatrix<AminoAcidCompound> getMatrix(String matrixName) {

		return matrices.get(matrixName);
		
	}
	
	public InputStream getInputStreamToAAindexFile(){
		 return SubstitutionMatrixHelper.class.getResourceAsStream(String.format("/AAINDEX.txt"));
	}

}
