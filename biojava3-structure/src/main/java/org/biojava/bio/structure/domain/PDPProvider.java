/**
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
 * Created on May 1, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.domain;

import java.util.SortedSet;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;

public interface PDPProvider {
	
	public SortedSet<String> getPDPDomainNamesForPDB(String pdbId);
	public Structure getDomain(String pdbDomainName, AtomCache cache);
	
}
