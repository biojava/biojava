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

import java.io.InputStream;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava3.protmod.io.ProteinModificationXmlReader;

/**
 * This class serves as a instance registry by maintaining 
 * a pool of ProteinModification instances.
 * 
 * A list of common protein modifications were preloaded
 * from an XML file.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ProteinModificationRegistry {
	private static Set<ProteinModification> registry = null;
	private static Map<String, ProteinModification> byId = null;
	private static Map<String, Set<ProteinModification>> byResidId = null;
	private static Map<String, Set<ProteinModification>> byPsimodId = null;
	private static Map<String, Set<ProteinModification>> byPdbccId = null;
	private static Map<String, Set<ProteinModification>> byKeyword = null;
	private static Map<Component, Set<ProteinModification>> byComponent = null;
	private static Map<ModificationCategory, Set<ProteinModification>> byCategory = null;
	private static Map<ModificationOccurrenceType, Set<ProteinModification>> byOccurrenceType = null;
	
	private static String DIR_XML_PTM_LIST = "ptm_list.xml";
	
	
	
	/**
	 * register common protein modifications from XML file.
	 */
	private static void registerCommonProteinModifications(InputStream inStream) {
		try {
			
			ProteinModificationXmlReader.registerProteinModificationFromXml(inStream);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Initialization the static variables and register common modifications.
	 */
	public static void init() {
		lazyInit();
		
	}
	
	/** Initialization the static variables and register common modifications.
	 * Allows external user to provide alternative ptm_list.xml file instead of the one contained in this jar file.
	 * 
	 * @param inStream InputStream to a XML file containing the list of PTMs (as in ptm_list.xml)
	 */
	 
	public static void init(InputStream inStream) {
		lazyInit(inStream);
	}
	
	
	
	/**
	 * Lazy Initialization the static variables and register common modifications.
	 * just opens the stream to ptm_list.xml and delegates to lazyInit(InputStream) for parsing. 
	 */
	private static synchronized void lazyInit() {
		if (registry==null) {	
			InputStream isXml = ProteinModification.class.getResourceAsStream(DIR_XML_PTM_LIST);
			lazyInit(isXml);
		}
	}
	
	
	/**
	 * Lazy Initialization the static variables and register common modifications. 
	 */
	private static synchronized void lazyInit(InputStream inStream) {
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
			registerCommonProteinModifications(inStream);					
		}
	}

	/**
	 * Register a new ProteinModification.
	 */
	public static void register(final ProteinModification modification) {
		if (modification==null) throw new IllegalArgumentException("modification == null!");
		
		lazyInit();

		String id = modification.getId();
		if (byId.containsKey(id)) {
			throw new IllegalArgumentException(id+" has already been registered.");
		}
		
		registry.add(modification);
		byId.put(id, modification);
		
		ModificationCategory cat = modification.getCategory();
		byCategory.get(cat).add(modification);
		
		ModificationOccurrenceType occType = modification.getOccurrenceType();
		byOccurrenceType.get(occType).add(modification);

		
		ModificationCondition condition = modification.getCondition();
		List<Component> comps = condition.getComponents();
		for (Component comp:comps) {
			Set<ProteinModification> mods = byComponent.get(comp);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byComponent.put(comp, mods);
			}
			mods.add(modification);
		}
		
		String pdbccId = modification.getPdbccId();
		if (pdbccId!=null) {
			Set<ProteinModification> mods = byPdbccId.get(pdbccId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byPdbccId.put(pdbccId, mods);
			}
			mods.add(modification);
		}
		
		String residId = modification.getResidId();
		if (residId!=null) {
			Set<ProteinModification> mods = byResidId.get(residId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byResidId.put(residId, mods);
			}
			mods.add(modification);
		}
		
		String psimodId = modification.getPsimodId();
		if (psimodId!=null) {			
			Set<ProteinModification> mods = byPsimodId.get(psimodId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byPsimodId.put(psimodId, mods);
			}
			mods.add(modification);
		}
		
		for (String keyword : modification.getKeywords()) {
			Set<ProteinModification> mods = byKeyword.get(keyword);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byKeyword.put(keyword, mods);
			}
			mods.add(modification);
		}
	}
	
	/**
	 * Remove a modification from registry.
	 * @param mod
	 */
	public static void unregister(ProteinModification modification) {
		if (modification==null) throw new IllegalArgumentException("modification == null!");
		
		registry.remove(modification);
		
		byId.remove(modification.getId());
		
		Set<ProteinModification> mods;
		
		mods = byResidId.get(modification.getResidId());
		if (mods!=null) mods.remove(modification);
		
		mods = byPsimodId.get(modification.getPsimodId());
		if (mods!=null) mods.remove(modification);
		
		mods = byPdbccId.get(modification.getPdbccId());
		if (mods!=null) mods.remove(modification);
		
		for (String keyword : modification.getKeywords()) {
			mods = byKeyword.get(keyword);
			if (mods!=null) mods.remove(modification);
		}
		
		ModificationCondition condition = modification.getCondition();
		List<Component> comps = condition.getComponents();
		for (Component comp : comps) {
			mods = byComponent.get(comp);
			if (mods!=null) mods.remove(modification);
		}
		
		byCategory.get(modification.getCategory()).remove(modification);
		byOccurrenceType.get(modification.getOccurrenceType()).remove(modification);
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
	
	
}
