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
 * Created on Nov 17, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.Collection;

/**
 * Uses ProteinModificationBuilder pattern to set optional attributes for a ProteinModification. 
 * For example, this allows you to use the following code:
 * <pre>
 * ProteinModification
 *     .register("0001", Modification.ATTACHMENT, 
 *     		ModificationOccurrenceType.NATURAL,
 *     		new ModificationCondition(Component.register("SER"),
 *     			Component.register("XYS"),"OG","C1"))
 *     .setResidId("AA0406")
 *     .setResidName("O-xylosyl-L-serine");
 * </pre>
 */
public final class ProteinModificationBuilder {
	private final ProteinModification current;
	
	/**
	 * Create a ProteinModificationBuilder for a ProteinModification. 
	 * @param current the ProteinModification to be modified.
	 */
	public ProteinModificationBuilder(final ProteinModification current) {
		if (current==null) {
			throw new IllegalArgumentException("Null argument.");
		}
		
		this.current = current;
	}
	
	/**
	 * 
	 * @return the ProteinModification under construction.
	 */
	public ProteinModification asModification() {
		return current;
	}
	
	/**
	 * Set the Protein Data Bank Chemical Component ID.
	 * @param pdbccId Protein Data Bank Chemical Component ID.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setPdbccId(final String pdbccId) {
		current.setPdbccId(pdbccId);
		return this;
	}
	
	/**
	 * Set the Protein Data Bank Chemical Component name.
	 * @param pdbccName Protein Data Bank Chemical Component name.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setPdbccName(final String pdbccName) {	
		current.setPdbccName(pdbccName);		
		return this;
	}
	
	/**
	 * Set the RESID ID.
	 * @param residId RESID ID.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setResidId(final String residId) {
		current.setResidId(residId);		
		return this;
	}
	
	/**
	 * Set the RESID name.
	 * @param residName RESID name.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setResidName(final String residName) {
		current.setResidName(residName);		
		return this;
	}
	
	/**
	 * Set the PSI-MOD ID.
	 * @param psimodId PSI-MOD ID.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setPsimodId(final String psimodId) {
		current.setPsimodId(psimodId);
		return this;
	}
	
	/**
	 * Set the PSI-MOD name.
	 * @param psimodName PSI-MOD name.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setPsimodName(final String psimodName) {
		current.setPsimodName(psimodName);		
		return this;
	}
	
	/**
	 * Set the systematic name.
	 * @param sysName systematic name.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setSystematicName(final String sysName) {	
		current.setSystematicName(sysName);		
		return this;
	}
	
	/**
	 * 
	 * @param description description of the modification.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setDescription(final String description) {
		current.setDescription(description);
		
		return this;
	}
	
	/**
	 * Add a keyword associate with the PTM.
	 * @param keyword a keyword.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 * @throws IllegalArgumentException if the keyword is null.
	 */
	public ProteinModificationBuilder addKeyword(String keyword) {
		if (keyword==null)
			throw new NullPointerException();
		current.addKeyword(keyword);
		return this;
	}
	
	public ProteinModificationBuilder addKeywords(Collection<String> keywords) {
		if (keywords==null) {
			throw new IllegalArgumentException("Keywords cannot be null.");
		}
		
		for (String keyword : keywords) {
			addKeyword(keyword);
		}
		
		return this;
	}
	
	/**
	 * Set the residue formula.
	 * @param formula residue formula.
	 * @return the same ProteinModificationBuilder object so you can chain setters.
	 */
	public ProteinModificationBuilder setFormula(final String formula) {
		current.setFormula(formula);		
		return this;
	}
}
