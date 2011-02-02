package org.biojava.bio.structure.io.mmcif;

import org.biojava.bio.structure.io.mmcif.model.ChemComp;

/** Interface that is implemented by all classes that can provide {@link ChemComp} definitions.
 * 
 * @author Andreas Prlic
 * @since 3.0
 */
public interface ChemCompProvider {

	/** Returns a new instance of a chemical component definition.
	 * 
	 * @param recordName the ID of the {@link ChemComp}
	 * @return a new {@link ChemComp} definition.
	 */
	ChemComp getChemComp(String recordName) ;

}
