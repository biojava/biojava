/**
 * 
 */
package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * CEHook which does no processing. This can be used as a base
 * class for other implementations.
 * 
 * @author Spencer Bliven
 *
 */
public class DefaultCEHook implements CEHook {

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.align.ce.CEHook#postprocess(org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.lang.Object, org.biojava.bio.structure.align.model.AFPChain)
	 */
	public AFPChain postprocess(Atom[] ca1m, Atom[] ca2m, Object param,
			AFPChain afpChain) {
		return afpChain;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.align.ce.CEHook#preprocess(org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.lang.Object)
	 */
	public Atom[][] preprocess(Atom[] ca1, Atom[] ca2, Object param) {
		return new Atom[][] {ca1,ca2};
	}

}
