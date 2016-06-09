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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A bean for identifying groups in GroupContactSets.
 * Used only within the contact package to be able to compare
 * contacts between chains of same entity/compound based on residue numbers
 * and independently from chain identifiers.
 *
 * @author duarte_j
 */
class ResidueIdentifier implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(ResidueIdentifier.class);

	private int seqResIndex;


	public ResidueIdentifier(Group g) {

		Chain c = g.getChain();
		if (c==null) {
			logger.warn("Chain is not available for group {}. Contact comparison will not work for this residue",g.toString());
			this.seqResIndex = -1;
		} else {
			EntityInfo comp = c.getEntityInfo();
			if (comp==null) {
				logger.warn("Compound is not available for group {}. Contact comparison will not work for this residue",g.toString());
				this.seqResIndex = -1;
			} else {
				this.seqResIndex = comp.getAlignedResIndex(g, c);
			}

		}
	}

	public int getSeqResIndex() {
		return seqResIndex;
	}

	public void setSeqResIndex(int seqResIndex) {
		this.seqResIndex = seqResIndex;
	}

	@Override
	public int hashCode() {
		return seqResIndex;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ResidueIdentifier other = (ResidueIdentifier) obj;

		return this.seqResIndex == other.seqResIndex;
	}

	@Override
	public String toString() {
		return String.valueOf(seqResIndex);
	}

}
