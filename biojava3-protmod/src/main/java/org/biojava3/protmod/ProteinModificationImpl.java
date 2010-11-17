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

import java.util.HashSet;
import java.util.Set;

/**
 * This class contains information about a specific protein
 * modification. 
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class ProteinModificationImpl 
implements ProteinModification {
	
	public ProteinModificationImpl() {
		
	}
	
	/**
	 * 
	 */
	public ProteinModificationImpl(final String id) {
		this.id = id;
	}
	
	/**
	 * 
	 * @param id
	 * @param cat
	 * @param occType
	 * @param condition
	 */
	public ProteinModificationImpl(final String id, final ModificationCategory cat,
			final ModificationOccurrenceType occType, final ModificationCondition condition) {
		this.id = id;
		this.category = cat;
		this.occurrenceType = occType;
		this.condition = condition;
		this.keywords = new HashSet<String>();
	}
	
	public ProteinModificationImpl(ProteinModification other) {
		this(other.getId(), other.getCategory(), other.getOccurrenceType(), other.getCondition());
		this.pdbccId = other.getPdbccId();
		this.pdbccName = other.getPdbccName();
		this.residId = other.getResidId();
		this.psimodId = other.getPsimodId();
		this.psimodName = other.getPsimodName();
		this.sysName = other.getSystematicName();
		this.formula = other.getFormula();
		this.description = other.getDescription();
		this.keywords = new HashSet<String>(other.getKeywords());
	}
	
	private String id;
	private String pdbccId = null;
	private String pdbccName = null;
	private String residId = null;
	private String residName = null;
	private String psimodId = null;
	private String psimodName = null;
	private String sysName = null;
	private String formula = null;
	private String description = null;
	private Set<String> keywords;
	
	private ModificationCondition condition = null;
	private ModificationCategory category = null;
	private ModificationOccurrenceType occurrenceType = null;

	public void setId(String id) {
		this.id = id;
	}

	public void setPdbccId(String pdbccId) {
		this.pdbccId = pdbccId;
	}

	public void setPdbccName(String pdbccName) {
		this.pdbccName = pdbccName;
	}

	public void setResidId(String residId) {
		this.residId = residId;
	}

	public void setResidName(String residName) {
		this.residName = residName;
	}

	public void setPsimodId(String psimodId) {
		this.psimodId = psimodId;
	}

	public void setPsimodName(String psimodName) {
		this.psimodName = psimodName;
	}

	public void setFormula(String formula) {
		this.formula = formula;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public void setSystematicName(String sysName) {
		this.sysName = sysName;
	}

	public void addKeyword(String keyword) {
		this.keywords.add(keyword);
	}

	public void setCondition(ModificationCondition condition) {
		this.condition = condition;
	}

	public void setCategory(ModificationCategory category) {
		this.category = category;
	}

	public void setOccurrenceType(ModificationOccurrenceType occurrenceType) {
		this.occurrenceType = occurrenceType;
	}

	/**
	 * 
	 * @return modification id.
	 */
	public String getId() {
		return id;
	}
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component ID.
	 */
	public String getPdbccId() {
		return pdbccId;
	}
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component name.
	 */
	public String getPdbccName() {
		return pdbccName;
	}
	
	/**
	 * 
	 * @return RESID ID.
	 */
	public String getResidId() {
		return residId;
	}
	
	/**
	 * 
	 * @return RESID name.
	 */
	public String getResidName() {
		return residName;
	}
	
	/**
	 * 
	 * @return PSI-MOD ID.
	 */
	public String getPsimodId() {
		return psimodId;
	}
	
	/**
	 * 
	 * @return PSI-MOD name.
	 */
	public String getPsimodName() {
		return psimodName;
	}
	
	/**
	 * 
	 * @return Systematic name.
	 */
	public String getSystematicName() {
		return sysName;
	}
	
	/**
	 * 
	 * @return Description.
	 */
	public String getDescription() {
		return description;
	}
	
	/**
	 * 
	 * @return a set of keywords.
	 */
	public Set<String> getKeywords() {
		return keywords;
	}
	
	/**
	 * 
	 * @return {@link ModificationCondition}
	 */
	public ModificationCondition getCondition() {
		return condition;
	}
	
	/**
	 * 
	 * @return formula of the modified residue.
	 */
	public String getFormula() {
		return formula;
	}
	
	/**
	 * 
	 * @return the modification category.
	 */
	public ModificationCategory getCategory() {
		return category;
	}
	
	/**
	 * 
	 * @return the modification occurrence type.
	 */
	public ModificationOccurrenceType getOccurrenceType() {
		return occurrenceType;
	}
	
	/**
	 * 
	 * @return informative description.
	 */
	@Override
	public String toString() {
		return printModification(this);
//		StringBuilder sb = new StringBuilder();
//		
//		sb.append("ID:"+getId());
//		sb.append("\nPDBCC ID:"+getPdbccId());
//		sb.append("\tPDBCC name:"+getPdbccName());
//		sb.append("\nRESID ID:"+getResidId());
//		sb.append("\tRESID name:"+getResidName());
//		sb.append("\nPSI-MOD ID:"+getPsimodId());
//		sb.append("\tPSI-MOD name:"+getPsimodName());
//		sb.append("\nDescription:"+getDescription());
//		sb.append("\nSystematic name:"+getSystematicName());
//		sb.append("\nCategory:"+getCategory().label());
//		sb.append("\nOccurrence type:"+getOccurrenceType().label());
//		sb.append("\nKeywords:"+getKeywords());
//		sb.append("\nCondition:"+getCondition());
//		return sb.toString();
	}
	
	private String printModification(ProteinModificationImpl mod) {
		StringBuilder sb = new StringBuilder();
						
		Set<String> keywords = mod.getKeywords();
		if (keywords!=null && !keywords.isEmpty()) {

			for (String keyword : keywords) {
				sb.append(keyword);
				sb.append(", ");
			}
			sb.delete(sb.length()-2,sb.length());
		}
		sb.append(" [ID:");
		sb.append(mod.getId());
		sb.append("] [");
		sb.append(mod.getCategory());
		sb.append("] ");
		
		
		String resid = mod.getResidId();
		if (resid != null) {
			sb.append("; ");
			sb.append("RESID:");
			sb.append(resid);
			String residname = mod.getResidName();
			if (residname != null) {
				sb.append(" (");
				sb.append(residname);
				sb.append(')');
			}
		}
		
		return sb.toString();
	}
	
	public int hashCode() {
		int ret = id.hashCode();
		ret = ret * 31 + category.hashCode();
		return ret;
	}
	
	public boolean equals(Object obj) {
		if (!(obj instanceof ProteinModification))
			return false;
		
		ProteinModification mod = (ProteinModification)obj;
		if (!id.equals(mod.getId()))
			return false;
		
		if (category != mod.getCategory())
			return false;
		
		return true;
	}
}
