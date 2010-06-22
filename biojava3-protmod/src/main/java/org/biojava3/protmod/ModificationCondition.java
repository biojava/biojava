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
 * Created on Jun 4, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.List;

/**
 * Conditions of a protein modification, e.g. components and atoms.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ModificationCondition {
	/**
	 * 
	 * @return the involved components.
	 */
	public List<Component> getComponents();
	
	/**
	 * 
	 * @return a list of all {Link ModificationLinkage}s.
	 */
	public List<ModificationLinkage> getLinkages();
}
