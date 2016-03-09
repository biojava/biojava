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

package org.biojava.nbio.ontology;




/**
 * Abstract implementation of term
 *
 * This provides basic change-forwarding functionality from
 *                the annotation and ontology properties.
 *
 * @author Thomas Down
 * @since 1.4
 */

public abstract class AbstractTerm  implements Term {

	protected String description;

	@Override
	public  void setDescription(String description){
		this.description = description;
	}
}
