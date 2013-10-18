/**
 * 
 */
package org.biojava3.core.sequence.io.template;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * @author mckeee1
 * 
 */
public interface GenbankHeaderFormatInterface<S extends Sequence<?>, C extends Compound> {
	/**
	 * 
	 * @param sequence
	 * @return
	 */
	public static final String UNKNOWN_DNA = "UNK";

	public String getHeader(S sequence);

}
