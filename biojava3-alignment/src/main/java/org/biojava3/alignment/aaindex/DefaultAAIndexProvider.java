package org.biojava3.alignment.aaindex;

import java.io.InputStream;

import java.util.Map;

import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;


/** The default provider for AAINDEX loads substituion matrices from the AAINDEX file in the resources directory
 * 
 * @author Andreas Prlic
 *
 */
public class DefaultAAIndexProvider implements AAIndexProvider {

	Map<String,SubstitutionMatrix<AminoAcidCompound>> matrices;
	
	public DefaultAAIndexProvider(){
		
		
		InputStream inStream = getInputStreamToAAindexFile();
		
		AAIndexFileParser parser = new AAIndexFileParser();
		
		try {
			parser.parse(inStream);
		} catch ( Exception e){
			e.printStackTrace();
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
