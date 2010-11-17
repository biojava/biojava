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
 * Created on May 27, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.Set;

/**
 * This interface defines information about a specific protein
 * modification. 
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ProteinModification {

	public void setId(String id);

	public void setPdbccId(String pdbccId);

	public void setPdbccName(String pdbccName);

	public void setResidId(String residId);

	public void setResidName(String residName);

	public void setPsimodId(String psimodId);

	public void setPsimodName(String psimodName);

	public void setFormula(String formula);

	public void setDescription(String description);

	public void setSystematicName(String sysName);

	public void addKeyword(String keyword);

	public void setCondition(ModificationCondition condition);

	public void setCategory(ModificationCategory category);

	public void setOccurrenceType(ModificationOccurrenceType occurrenceType);

	/**
	 * 
	 * @return modification id.
	 */
	public String getId();
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component ID.
	 */
	public String getPdbccId();
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component name.
	 */
	public String getPdbccName();
	
	/**
	 * 
	 * @return RESID ID.
	 */
	public String getResidId();
	
	/**
	 * 
	 * @return RESID name.
	 */
	public String getResidName();
	
	/**
	 * 
	 * @return PSI-MOD ID.
	 */
	public String getPsimodId();
	
	/**
	 * 
	 * @return PSI-MOD name.
	 */
	public String getPsimodName();
	
	/**
	 * 
	 * @return Systematic name.
	 */
	public String getSystematicName();
	
	/**
	 * 
	 * @return Description.
	 */
	public String getDescription();
	
	/**
	 * 
	 * @return a set of keywords.
	 */
	public Set<String> getKeywords();
	
	/**
	 * 
	 * @return {@link ModificationCondition}
	 */
	public ModificationCondition getCondition();
	
	/**
	 * 
	 * @return formula of the modified residue.
	 */
	public String getFormula();
	
	/**
	 * 
	 * @return the modification category.
	 */
	public ModificationCategory getCategory();
	
	/**
	 * 
	 * @return the modification occurrence type.
	 */
	public ModificationOccurrenceType getOccurrenceType();
	
}
