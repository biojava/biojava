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

import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * This class contains information about a specific protein
 * modification. 
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationImpl 
implements ProteinModification , Comparable<ProteinModification> {
	
	private final String id;	
	private final ModificationCondition condition;
	private final ModificationCategory category;
	private final ModificationOccurrenceType occurrenceType;
	
	private final String pdbccId;
	private final String pdbccName;
	private final String residId;
	private final String residName;
	private final String psimodId;
	private final String psimodName;
	private final String sysName;
	private final String formula;
	
	private final Set<String> keywords;

	public String getId() {
		return id;
	}
	
	public String getPdbccId() {
		return pdbccId;
	}
	
	public String getPdbccName() {
		return pdbccName;
	}
	
	public String getResidId() {
		return residId;
	}
	
	public String getResidName() {
		return residName;
	}
	
	public String getPsimodId() {
		return psimodId;
	}
	
	public String getPsimodName() {
		return psimodName;
	}
	
	public String getSystematicName() {
		return sysName;
	}
	
	public String getDescription() {
		return toString();
		//return description;
	}
	
	public Set<String> getKeywords() {
		return keywords;
	}
	
	public ModificationCondition getCondition() {
		return condition;
	}
	
	public String getFormula() {
		return formula;
	}
	
	public ModificationCategory getCategory() {
		return category;
	}
	
	public ModificationOccurrenceType getOccurrenceType() {
		return occurrenceType;
	}
	
	@Override
	public String toString() {
		return "ProteinModificationImpl [id=" + id + ", condition=" + condition
				+ ", category=" + category + ", occurrenceType="
				+ occurrenceType + ", pdbccId=" + pdbccId + ", pdbccName="
				+ pdbccName + ", residId=" + residId + ", residName="
				+ residName + ", psimodId=" + psimodId + ", psimodName="
				+ psimodName + ", sysName=" + sysName + ", formula=" + formula
				+ ", keywords=" + keywords + "]";
	}
	
	private String printModification(ProteinModificationImpl mod) {
		StringBuilder sb = new StringBuilder();
						
		String name = getBestPossibleName(mod);
		boolean hasName = true;
		if (  name.equals(""))
			hasName = false;
		sb.append(name);
		
		Set<String> keywords = mod.getKeywords();
		if (keywords!=null && !keywords.isEmpty()) {
			if ( hasName)
				sb.append(" (");
			for (String keyword : keywords) {
				
				sb.append(keyword);
				sb.append(", ");
			}
			sb.delete(sb.length()-2,sb.length());
		}
		if ( hasName)
			sb.append(")");
		return sb.toString();
	}
	
	private String getBestPossibleName(ProteinModificationImpl mod) {
		
		//System.out.println(mod.getResidName() + " : " + mod.getPsimodName() + " : " + mod.getPdbccName() + " : " + mod.getSystematicName());
		
		// first: get resid
		String resid = mod.getResidId();
		if (resid != null) {
			String residname = mod.getResidName();
			if (residname != null) {
				return residname;
			}
		}
		
		// 2nd: PSI-MOD
		
		String name = mod.getPsimodName();
		if ( name != null) {
			//System.out.println("PSI_MOD name:" + name);
			return name;
		}
		
		// 3rd PDB-CC
		
		String pdbcc = mod.getPdbccName();
		if ( pdbcc != null ) {
			//System.out.println("PDBCC name: " + pdbcc);
			return pdbcc;
		}
		
		
		// no public name know, use the systematic name
		
		String systematic = mod.getSystematicName();
		if ( systematic != null) {
			//System.out.println("SYSTEMATIC NAME: " + mod.getSystematicName());
			return systematic;
		}
		
		
		return "";
		
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
	
	

	
	/**
	 * Uses Builder pattern to build a ProteinModification.
	 */
	public static class Builder {
		private final String id;
		
		private ModificationCondition condition;
		private ModificationCategory category;
		private ModificationOccurrenceType occurrenceType;
		
		private String pdbccId = null;
		private String pdbccName = null;
		private String residId = null;
		private String residName = null;
		private String psimodId = null;
		private String psimodName = null;
		private String sysName = null;
		private String formula = null;
		
		private Set<String> keywords = new LinkedHashSet<String>();
		
		/**
		 * 
		 * @param id
		 * @param cat
		 * @param occType
		 * @param condition
		 */
		public Builder(final String id, final ModificationCategory cat,
				final ModificationOccurrenceType occType, 
				final ModificationCondition condition) {
			if ( id == null) throw new IllegalArgumentException("id == null!");
			if ( cat == null) throw new IllegalArgumentException("cat == null!");
			if ( occType == null) throw new IllegalArgumentException("occType == null!");
			if ( condition == null) throw new IllegalArgumentException("condition == null!");
			
			this.id = id;
			this.category = cat;
			this.occurrenceType = occType;
			this.condition = condition;
		}
		
		/**
		 * Create a Builder from an existing ProteinModification. 
		 * @param copyFrom the ProteinModification to be copied from.
		 */
		public Builder(final ProteinModification copyFrom) {
			this(copyFrom.getId(), copyFrom.getCategory(), copyFrom.getOccurrenceType(), copyFrom.getCondition());
			this.pdbccId = copyFrom.getPdbccId();
			this.pdbccName = copyFrom.getPdbccName();
			this.residId = copyFrom.getResidId();
			this.residName = copyFrom.getResidName();
			this.psimodId = copyFrom.getPsimodId();
			this.psimodName = copyFrom.getPsimodName();
			this.sysName = copyFrom.getSystematicName();
			this.formula = copyFrom.getFormula();
			
			this.keywords = new LinkedHashSet<String>(copyFrom.getKeywords());
		}
		
		public Builder setCategory(final ModificationCategory cat) {
			if (cat == null) throw new IllegalArgumentException("cat == null!");
			this.category = cat;
			return this;
		}
		
		public Builder setOccurrenceType(final ModificationOccurrenceType occType) {
			if (occType == null) throw new IllegalArgumentException("occType == null!");
			this.occurrenceType =occType;
			return this;
		}
		
		public Builder setCondition(final ModificationCondition condition) {
			if (condition == null) throw new IllegalArgumentException("condition == null!");
			this.condition = condition;
			return this;
		}
		
		/**
		 * Set the Protein Data Bank Chemical Component ID.
		 * @param pdbccId Protein Data Bank Chemical Component ID.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setPdbccId(final String pdbccId) {
			this.pdbccId = pdbccId;
			return this;
		}
		
		/**
		 * Set the Protein Data Bank Chemical Component name.
		 * @param pdbccName Protein Data Bank Chemical Component name.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setPdbccName(final String pdbccName) {	
			this.pdbccName = pdbccName;		
			return this;
		}
		
		/**
		 * Set the RESID ID.
		 * @param residId RESID ID.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setResidId(final String residId) {
			this.residId = residId;		
			return this;
		}
		
		/**
		 * Set the RESID name.
		 * @param residName RESID name.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setResidName(final String residName) {
			this.residName = residName;		
			return this;
		}
		
		/**
		 * Set the PSI-MOD ID.
		 * @param psimodId PSI-MOD ID.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setPsimodId(final String psimodId) {
			this.psimodId = psimodId;
			return this;
		}
		
		/**
		 * Set the PSI-MOD name.
		 * @param psimodName PSI-MOD name.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setPsimodName(final String psimodName) {
			this.psimodName = psimodName;		
			return this;
		}
		
		/**
		 * Set the systematic name.
		 * @param sysName systematic name.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setSystematicName(final String sysName) {	
			this.sysName = sysName;		
			return this;
		}
		
		/**
		 * 
		 * @param description description of the modification.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setDescription(final String description) {
			// description is created on the fly in getDescription
			return this;
		}
		
		/**
		 * Add a keyword associate with the PTM.
		 * @param keyword a keyword.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if the keyword is null.
		 */
		public Builder addKeyword(String keyword) {
			if (keyword == null) throw new IllegalArgumentException("Keyword cannot be null.");
			keywords.add(keyword);
			return this;
		}
		
		public Builder addKeywords(Collection<String> keywords) {
			if (keywords==null)	throw new IllegalArgumentException("Keywords cannot be null.");
			
			for (String keyword : keywords) {
				addKeyword(keyword);
			}
			
			return this;
		}
		
		/**
		 * Set the residue formula.
		 * @param formula residue formula.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder setFormula(final String formula) {
			this.formula = formula;		
			return this;
		}
		
		/**
		 * 
		 * @return build ProteinModification.
		 */
		public ProteinModificationImpl build() {
			return new ProteinModificationImpl(this);
		}
	}
	
	/**
	 * 
	 */
	private ProteinModificationImpl(Builder builder) {
		this.id = builder.id;
		this.category = builder.category;
		this.occurrenceType = builder.occurrenceType;
		this.condition = builder.condition;
		this.pdbccId = builder.pdbccId;
		this.pdbccName = builder.pdbccName;
		this.residId = builder.residId;
		this.residName = builder.residName;
		this.psimodId = builder.psimodId;
		this.psimodName = builder.psimodName;
		this.sysName = builder.sysName;
		this.formula = builder.formula;
		
		this.keywords = new LinkedHashSet<String>(builder.keywords);
	}

	@Override
	public int compareTo(ProteinModification arg0) {
		if ( this.equals(arg0))
			return 0;
		
		return this.toString().compareTo(arg0.toString());
	}
}
