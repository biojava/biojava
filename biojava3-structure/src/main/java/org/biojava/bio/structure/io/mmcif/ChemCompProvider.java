package org.biojava.bio.structure.io.mmcif;

import org.biojava.bio.structure.io.mmcif.model.ChemComp;

/** Interface that is implememted by all classes that can instantiate ChemComp groups
 * 
 * @author Andreas Prlic
 *
 */
public interface ChemCompProvider {

	ChemComp getChemComp(String recordName) ;

}
