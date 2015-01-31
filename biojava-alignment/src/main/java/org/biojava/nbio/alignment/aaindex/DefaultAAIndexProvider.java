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
package org.biojava.nbio.alignment.aaindex;

import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;


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
