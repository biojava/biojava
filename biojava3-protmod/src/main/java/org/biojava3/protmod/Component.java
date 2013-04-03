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

import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * contains information about a certain Component.
 * The Component class uses the extensible enum pattern.
 * You can't instantiate Component directly, instead 
 * you have to use one of the {@link register} and {@link of} methods.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class Component {
	private final Set<String> pdbccIds;
	private final boolean isNTerminal;
	private final boolean isCTerminal;
	
	private static Set<Component> components = null;
	private static Map<Set<String>, Component> nonTerminalComps = null;
	private static Map<Set<String>, Component> nTerminalAminoAcids = null;
	private static Map<Set<String>, Component> cTerminalAminoAcids = null;
	
	/**
	 * Lazy initialization of the static variables.
	 */
	private static void lazyInit() {
		if (components==null) {
			components = new HashSet<Component>();
			nonTerminalComps = new HashMap<Set<String>, Component>();
			nTerminalAminoAcids = new HashMap<Set<String>, Component>();
			cTerminalAminoAcids = new HashMap<Set<String>, Component>();
		}
	}
	
	/**
	 * Create a ComponentImpl. 
	 * @param pdbccIds a set of possible Protein Data Bank ID. Cannot be null or empty.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or terminal condition is indicated for non-amino-acid component,
	 *  or both N-terminal and C-terminal are true.
	 */
	private Component(final Set<String> pdbccIds, 
			final boolean isNTerminal, final boolean isCTerminal) {
		if (pdbccIds==null || pdbccIds.isEmpty()) {
			throw new IllegalArgumentException("pdbccId or type cannot be null.");
		}
		
		if (isNTerminal&&isCTerminal) {
			throw new IllegalArgumentException("An amino acid can be specified at" +
					"N-terminal or C-terminal but not both."); //TODO: is this true?
		}
		
		this.pdbccIds = pdbccIds;
		this.isNTerminal = isNTerminal;
		this.isCTerminal = isCTerminal;
	}
	
	/**
	 * 
	 * @return Protein Data Bank ID.
	 */
	public Set<String> getPdbccIds() {
		return pdbccIds;
	}
	
	/**
	 * 
	 * @return true if occurring on N terminal; false, otherwise.
	 */
	public boolean isNTerminal() {
		return isNTerminal;
	}

	/**
	 * 
	 * @return true if occurring on C terminal; false, other wise.
	 */
	public boolean isCTerminal() {
		return isCTerminal;
	}
	
	/**
	 * Get a Component that does not have to occur at terminals. If the 
	 * corresponding component has already been registered, return that one.
	 * @param pdbccIds possible Protein Data Bank ID.
	 * @return a component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type.
	 */
	public static Component of(final String pdbccId) {
		return of(pdbccId, false, false);
	}
	
	/**
	 * Get or create a Component.
	 * @param pdbccId Protein Data Bank ID.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @return a component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type,
	 *  or terminal condition is indicated for non-amino-acid component,
	 *  or both N-terminal and C-terminal are true.
	 */
	public static Component of(final String pdbccId, 
			final boolean isNTerminal, final boolean isCTerminal) {
		return of (Collections.singleton(pdbccId), isNTerminal, isCTerminal);
	}
	
	/**
	 * Get a Component that does not have to occur at terminals. If the 
	 * corresponding component has already been registered, return that one.
	 * @param pdbccIds a set of possible Protein Data Bank ID.
	 * @return a component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type.
	 */
	public static Component of(final Set<String> pdbccIds) {
		return of(pdbccIds, false, false);
	}
	
	/**
	 * Get or create a Component.
	 * @param pdbccIds a set of possible Protein Data Bank ID.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @return a component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type,
	 *  or terminal condition is indicated for non-amino-acid component,
	 *  or both N-terminal and C-terminal are true.
	 */
	public static Component of(final Set<String> pdbccIds, 
			final boolean isNTerminal, final boolean isCTerminal) {
		if (isNTerminal && isCTerminal) {
			throw new IllegalArgumentException("An amino acid can be at" +
			"N-terminal or C-terminal but not both."); //TODO: is this true?
		}
		
		lazyInit();
		
		Component comp;
		if (isNTerminal) {
			comp = nTerminalAminoAcids.get(pdbccIds);
			if (comp == null) {
				comp = new Component(pdbccIds, isNTerminal, isCTerminal);
				nTerminalAminoAcids.put(pdbccIds, comp);
			}
		} else if (isCTerminal) {
			comp = cTerminalAminoAcids.get(pdbccIds);
			if (comp == null) {
				comp = new Component(pdbccIds, isNTerminal, isCTerminal);
				cTerminalAminoAcids.put(pdbccIds, comp);
			}
		} else {
			comp = nonTerminalComps.get(pdbccIds);
			if (comp == null) {
				comp = new Component(pdbccIds, isNTerminal, isCTerminal);
				nonTerminalComps.put(pdbccIds, comp);
			}
		}
		
		components.add(comp);
		return comp;
	}
	
	public static Set<Component> allComponents() {
		return Collections.unmodifiableSet(components);
	}
	
	/**
	 * 
	 * @return informative description.
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getPdbccIds());
		if (isCTerminal()) {
			sb.append("(C)");
		} else if (isNTerminal()) {
			sb.append("(N)");
		}
		return sb.toString();
	}
}
