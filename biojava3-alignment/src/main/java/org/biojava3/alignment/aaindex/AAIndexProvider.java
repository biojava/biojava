package org.biojava3.alignment.aaindex;

import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.compound.AminoAcidCompound;


public interface AAIndexProvider {
	
	public SubstitutionMatrix<AminoAcidCompound> getMatrix(String matrixName);
}
