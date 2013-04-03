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

import java.io.InputStream;

import java.util.Collection;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava3.protmod.io.ProteinModificationXmlReader;

/**
 * This class contains information about a specific protein
 * modification. By applying extensible enum pattern, this
 * class also serves as a instance registry by maintaining 
 * a pool of ProteinModification instances. You can't 
 * instantiate ProteinModifications directly, instead you 
 * have to use {@link register} to create and register a
 * ProteinModification, and use one of the getBy... methods
 * to retrieve a ProteinModification.
 * 
 * A list of common protein modifications were preloaded
 * from an XML file.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class ProteinModification {
	
	/**
	 * Constructor is private, so that we don't
	 * get any stand-alone ProteinModifications. 
	 * ProteinModifications should be obtained from 
	 * getBy... methods. Information about
	 * DataSources can be added with {@link register}.
	 */
	private ProteinModification() {}
	
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
	private Set<String> keywords = new HashSet<String>();
	
	private ModificationCondition condition = null;
	
	private ModificationCategory category;
	private ModificationOccurrenceType occurrenceType;
	
	private static Set<ProteinModification> registry = null;
	private static Map<String, ProteinModification> byId = null;
	private static Map<String, Set<ProteinModification>> byResidId = null;
	private static Map<String, Set<ProteinModification>> byPsimodId = null;
	private static Map<String, Set<ProteinModification>> byPdbccId = null;
	private static Map<String, Set<ProteinModification>> byKeyword = null;
	private static Map<Component, Set<ProteinModification>> byComponent = null;
	private static Map<ModificationCategory, Set<ProteinModification>> byCategory = null;
	private static Map<ModificationOccurrenceType, Set<ProteinModification>> byOccurrenceType = null;
	
	/**
	 * Initialization the static variables and register common modifications.
	 */
	public static void init() {
		lazyInit();
	}
	
	/**
	 * Lazy Initialization the static variables and register common modifications. 
	 */
	private static void lazyInit() {
		if (registry==null) {	
			registry = 	new HashSet<ProteinModification>();
			byId = new HashMap<String, ProteinModification>();
			byResidId = new HashMap<String, Set<ProteinModification>>();
			byPsimodId = new HashMap<String, Set<ProteinModification>>();
			byPdbccId = new HashMap<String, Set<ProteinModification>>();
			byKeyword = new HashMap<String, Set<ProteinModification>>();
			byComponent = new HashMap<Component, Set<ProteinModification>>();
			byCategory = new EnumMap<ModificationCategory, Set<ProteinModification>>(
					ModificationCategory.class);
			for (ModificationCategory cat:ModificationCategory.values()) {
				byCategory.put(cat, new HashSet<ProteinModification>());
			}
			byOccurrenceType = new EnumMap<ModificationOccurrenceType, Set<ProteinModification>>(
					ModificationOccurrenceType.class);
			for (ModificationOccurrenceType occ:ModificationOccurrenceType.values()) {
				byOccurrenceType.put(occ, new HashSet<ProteinModification>());
			}
			registerCommonProteinModifications();
		}
	}
	
	private static String DIR_XML_PTM_LIST = "ptm_list.xml";
	
	/**
	 * register common protein modifications from XML file.
	 */
	private static void registerCommonProteinModifications() {
		try {
			InputStream isXml = ProteinModification.class.getResourceAsStream(DIR_XML_PTM_LIST);
			ProteinModificationXmlReader.registerProteinModificationFromXml(isXml);
		} catch (Exception e) {
			e.printStackTrace();
		}
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
	 * Uses builder pattern to set optional attributes for a ProteinModification. 
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
	public static final class Builder {
		private final ProteinModification current;
		
		/**
		 * Create a Builder for a ProteinModification. 
		 * This constructor should only be called by the register method.
		 * @param current the ProteinModification to be modified.
		 */
		private Builder(final ProteinModification current) {
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
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if pdbccId has been set.
		 */
		public Builder setPdbccId(final String pdbccId) {			
			if (current.pdbccId!=null) {
				throw new IllegalArgumentException("PDBCC ID has been set.");
			}
			
			current.pdbccId = pdbccId;
			
			Set<ProteinModification> mods = byPdbccId.get(pdbccId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byPdbccId.put(pdbccId, mods);
			}
			mods.add(current);
			
			return this;
		}
		
		/**
		 * Set the Protein Data Bank Chemical Component name.
		 * @param pdbccName Protein Data Bank Chemical Component name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if pdbccName has been set.
		 */
		public Builder setPdbccName(final String pdbccName) {			
			if (current.pdbccName!=null) {
				throw new IllegalArgumentException("PDBCC name has been set.");
			}
			
			current.pdbccName = pdbccName;
			
			return this;
		}
		
		/**
		 * Set the RESID ID.
		 * @param residId RESID ID.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if residId is null or
		 *  it has been set.
		 * @throws IllegalArgumentException if residIdhas been set.
		 */
		public Builder setResidId(final String residId) {
			if (current.residId!=null) {
				throw new IllegalArgumentException("RESID ID has been set.");
			}
			
			current.residId = residId;
			
			Set<ProteinModification> mods = byResidId.get(residId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byResidId.put(residId, mods);
			}
			mods.add(current);
			
			return this;
		}
		
		/**
		 * Set the RESID name.
		 * @param residName RESID name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if residId has been set.
		 */
		public Builder setResidName(final String residName) {
			if (current.residName!=null) {
				throw new IllegalArgumentException("RESID name has been set.");
			}
			
			current.residName = residName;
			
			return this;
		}
		
		/**
		 * Set the PSI-MOD ID.
		 * @param psimodId PSI-MOD ID.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if psimodId has been set.
		 */
		public Builder setPsimodId(final String psimodId) {
			if (current.psimodId!=null) {
				throw new IllegalArgumentException("PSI-MOD ID has been set.");
			}
			
			current.psimodId = psimodId;
			
			Set<ProteinModification> mods = byPsimodId.get(psimodId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byPsimodId.put(psimodId, mods);
			}
			mods.add(current);
			
			return this;
		}
		
		/**
		 * Set the PSI-MOD name.
		 * @param psimodName PSI-MOD name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if psimodName has been set.
		 */
		public Builder setPsimodName(final String psimodName) {
			if (current.psimodName!=null) {
				throw new IllegalArgumentException("PSI-MOD name has been set.");
			}
			
			current.psimodName = psimodName;
			
			return this;
		}
		
		/**
		 * Set the systematic name.
		 * @param sysName systematic name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if sysName has been set.
		 */
		public Builder setSystematicName(final String sysName) {			
			if (current.sysName!=null) {
				throw new IllegalArgumentException("Systematic name has been set.");
			}
			
			current.sysName = sysName;
			
			return this;
		}
		
		/**
		 * 
		 * @param description description of the modification.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if description has been set.
		 */
		public Builder setDescription(final String description) {
			if (current.description!=null) {
				throw new IllegalArgumentException("Description has been set.");
			}
			
			current.description = description;
			
			return this;
		}
		
		/**
		 * Add a keyword associate with the PTM.
		 * @param keyword a keyword.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if the keyword is null.
		 */
		public Builder addKeyword(String keyword) {
			if (keyword==null) {
				throw new IllegalArgumentException("Keyword cannot be null.");
			}
			
			current.keywords.add(keyword);
			
			Set<ProteinModification> mods = byKeyword.get(keyword);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byKeyword.put(keyword, mods);
			}
			mods.add(current);
			
			return this;
		}
		
		public Builder addKeywords(Collection<String> keywords) {
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
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if formula has been set.
		 */
		public Builder setFormula(final String formula) {
			if (current.formula!=null) {
				throw new IllegalArgumentException("Formula has been set.");
			}
			
			current.formula = formula;
			
			return this;
		}
	}

	/**
	 * Register a new ProteinModification with (optional) detailed information.
	 * After registration, you could set the optional parameters as the following code:
	 * <pre>
	 * ProteinModification
	 *     .register("0001", Modification.ATTACHMENT, 
	 *     		ModificationOccurrenceType.NATURAL,
	 *     		new ModificationCondition(Component.register("SER"),
	 *     			Component.register("XYS"),"OG","C1"))
	 *     .setResidId("AA0406")
	 *     .setResidName("O-xylosyl-L-serine");
	 * </pre>
	 *
	 * @param id modification id.
	 * @param cat modification category.
	 * @param occType occurrence type.
	 * @param condition {@link ModificationCondition}.
	 * @return Builder that can be used for adding detailed information.
	 * @throws IllegalArgumentException if null argument(s), or components is null
	 *  or empty in condition.
	 */
	public static Builder register(final String id, final ModificationCategory cat,
			final ModificationOccurrenceType occType, final ModificationCondition condition) {
		if ( id == null) throw new IllegalArgumentException("id == null!");
		if ( cat == null) throw new IllegalArgumentException("cat == null!");
		if ( occType == null) throw new IllegalArgumentException("occType == null!");
		if ( condition == null) throw new IllegalArgumentException("condition == null!");
		if (id==null || cat==null || occType==null || condition==null) {
			throw new IllegalArgumentException("Null argument(s)!");
		}
		
		List<Component> comps = condition.getComponents();
		if (comps==null || comps.isEmpty()) {
			throw new IllegalArgumentException("At least one component for" +
				" a modification.");
		}
		
		lazyInit();
		
		if (byId.containsKey(id)) {
			throw new IllegalArgumentException(id+" has already been registered.");
		}
		
		ProteinModification current = new ProteinModification();
		current.id = id;
		current.category = cat;
		current.occurrenceType = occType;
		current.condition = condition;		
		
		registry.add(current);
		byId.put(id, current);
		byCategory.get(cat).add(current);
		byOccurrenceType.get(occType).add(current);
		
		for (Component comp:comps) {
			Set<ProteinModification> mods = byComponent.get(comp);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byComponent.put(comp, mods);
			}
			mods.add(current);
		}
		
		return new Builder(current);
	}
	
	/**
	 * 
	 * @param id modification ID.
	 * @return ProteinModification that has the corresponding ID.
	 */
	public static ProteinModification getById(final String id) {	
		lazyInit();
		return byId.get(id);
	}

	/**
	 * 
	 * @param residId RESID ID.
	 * @return a set of ProteinModifications that have the RESID ID.
	 */
	public static Set<ProteinModification> getByResidId(final String residId) {
		lazyInit();
		return byResidId.get(residId);
	}
	/**
	 * 
	 * @param psimodId PSI-MOD ID.
	 * @return a set of ProteinModifications that have the PSI-MOD ID.
	 */
	public static Set<ProteinModification> getByPsimodId(final String psimodId) {
		lazyInit();
		return byPsimodId.get(psimodId);
	}
	
	/**
	 * 
	 * @param pdbccId Protein Data Bank Chemical Component ID.
	 * @return a set of ProteinModifications that have the PDBCC ID.
	 */
	public static Set<ProteinModification> getByPdbccId(final String pdbccId) {
		lazyInit();
		return byPdbccId.get(pdbccId);
	}
	
	/**
	 * 
	 * @param keyword a keyword.
	 * @return a set of ProteinModifications that have the keyword.
	 */
	public static Set<ProteinModification> getByKeyword(final String keyword) {
		lazyInit();
		return byKeyword.get(keyword);
	}
	
	/**
	 * Get ProteinModifications that involves one or more components.
	 * @param comp1 a {@link Component}.
	 * @param comps other {@link Component}s.
	 * @return a set of ProteinModifications that involves all the components.
	 */
	public static Set<ProteinModification> getByComponent(final Component comp1,
			final Component... comps) {
		lazyInit();
		Set<ProteinModification> mods = byComponent.get(comp1);
		if (mods==null) {
			return Collections.emptySet();
		}
		
		if (comps.length==0) {
			return Collections.unmodifiableSet(mods);
		} else {
			Set<ProteinModification> ret = new HashSet<ProteinModification>(mods);
			for (Component comp:comps) {
				mods = byComponent.get(comp);
				if (mods==null) {
					return Collections.emptySet();
				} else {
					ret.retainAll(mods);
				}
			}
			
			return ret;
		}
	}
	
	/**
	 * 
	 * @return set of all registered ProteinModifications.
	 */
	public static Set<ProteinModification> allModifications() {
		lazyInit();
		return Collections.unmodifiableSet(registry);
	}
	
	/**
	 * 
	 * @param cat {@link ModificationCategory}.
	 * @return set of registered ProteinModifications in a particular category.
	 */
	public static Set<ProteinModification> getByCategory(final ModificationCategory cat) {
		lazyInit();
		Set<ProteinModification> ret = byCategory.get(cat);
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @param occ {@link ModificationOccurrenceType}.
	 * @return set of registered ProteinModifications of a particular occurrence type.
	 */
	public static Set<ProteinModification> getByOccurrenceType(final ModificationOccurrenceType occ) {
		lazyInit();
		Set<ProteinModification> ret = byOccurrenceType.get(occ);
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of IDs of all registered ProteinModifications.
	 */
	public static Set<String> allIds() {
		lazyInit();
		Set<String> ret = byId.keySet();
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of PDBCC IDs of all registered ProteinModifications.
	 */
	public static Set<String> allPdbccIds() {
		lazyInit();
		Set<String> ret = byPdbccId.keySet();
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of RESID IDs of all registered ProteinModifications.
	 */
	public static Set<String> allResidIds() {
		lazyInit();
		Set<String> ret = byResidId.keySet();
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of PSI-MOD IDs of all registered ProteinModifications.
	 */
	public static Set<String> allPsimodIds() {
		lazyInit();
		Set<String> ret = byPsimodId.keySet();
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of components involved in all registered ProteinModifications.
	 */
	public static Set<Component> allComponents() {
		lazyInit();
		Set<Component> ret = byComponent.keySet();
		return Collections.unmodifiableSet(ret);
	}
	
	/**
	 * 
	 * @return set of keywords of all registered ProteinModifications.
	 */
	public static Set<String> allKeywords() {
		lazyInit();
		Set<String> ret = byKeyword.keySet();
		return Collections.unmodifiableSet(ret);
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
	
	private String printModification(ProteinModification mod) {
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
}
