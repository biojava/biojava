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

public class ComponentImpl implements Component {
	private final String pdbccId;
	private final boolean isNTerminal;
	private final boolean isCTerminal;
	private final ComponentType type;
	
	/**
	 * Create a ComponentImpl. Use this constructor if the component is
	 * not an amino acid or it does not occur at either terminal.
	 * @param pdbccId Protein Data Bank ID. Cannot be null.
	 * @param type {@link ComponentType}. Cannot be null.
	 * @throws IllegalArgumentException if pdbccId or type is null.
	 */
	public ComponentImpl(final String pdbccId, final ComponentType type) {
		this(pdbccId, type, false, false);
	}
	
	/**
	 * Create a ComponentImpl. Use this constructor if the component is
	 * an amino acid and occurs at either terminal.
	 * @param pdbccId Protein Data Bank ID. Cannot be null.
	 * @param type {@link ComponentType}. Cannot be null.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or terminal condition is indicated for non-amino-acid component,
	 *  or both N-terminal and C-terminal are true.
	 */
	public ComponentImpl(final String pdbccId, final ComponentType type,
			final boolean isNTerminal, final boolean isCTerminal) {
		if (pdbccId==null || type==null) {
			throw new IllegalArgumentException("pdbccId or type cannot be null.");
		}
		
		if ((isNTerminal||isCTerminal)&&(type!=ComponentType.AMINOACID)) {
			throw new IllegalArgumentException("Only amino acid can be specified" +
					" as N/C-terminal.");
		}
		
		if (isNTerminal&&isCTerminal) {
			throw new IllegalArgumentException("An amino acid can be specified at" +
					"N-terminal or C-terminal but not both."); //TODO: is this true?
		}
		
		this.pdbccId = pdbccId;
		this.type = type;
		this.isNTerminal = isNTerminal;
		this.isCTerminal = isCTerminal;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public String getPdbccId() {
		return pdbccId;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public boolean isNTerminal() {
		return isNTerminal;
	}

	/**
	 * {@inheritDoc}
	 */
	public boolean isCTerminal() {
		return isCTerminal;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public ComponentType getType() {
		return type;
	}
	
	/**
	 * 
	 * @return 
	 */
	public String toString() {
		return pdbccId;
	}
}
