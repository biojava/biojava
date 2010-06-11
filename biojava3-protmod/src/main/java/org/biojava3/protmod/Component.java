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

import java.util.HashMap;
import java.util.Map;

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
	private final String pdbccId;
	private final boolean isNTerminal;
	private final boolean isCTerminal;
	private final ComponentType type;
	
	private static Map<String, Component> nonTerminalComps = null;
	private static Map<String, Component> nTerminalAminoAcids = null;
	private static Map<String, Component> cTerminalAminoAcids = null;
	
	/**
	 * Lazy initialization of the static variables.
	 */
	private static void lazyInit() {
		if (nonTerminalComps==null) {
			nonTerminalComps = new HashMap<String, Component>();
			nTerminalAminoAcids = new HashMap<String, Component>();
			cTerminalAminoAcids = new HashMap<String, Component>();
		}
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
	private Component(final String pdbccId, final ComponentType type,
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
	 * 
	 * @return Protein Data Bank ID.
	 */
	public String getPdbccId() {
		return pdbccId;
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
	 * 
	 * @return the component type.
	 */
	public ComponentType getType() {
		return type;
	}
	
	/**
	 * Register a non-terminal Component. If the corresponding component 
	 * has already been registered, return that one.
	 * @param pdbccId Protein Data Bank ID. Cannot be null.
	 * @param type {@link ComponentType}. Cannot be null.
	 * @return the registered component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type.
	 */
	public static Component register(final String pdbccId, final ComponentType type) {
		return register(pdbccId, type, false, false);
	}
	
	/**
	 * Register a Component. If the corresponding component has already been 
	 * registered, return that one.
	 * @param pdbccId Protein Data Bank ID. Cannot be null.
	 * @param type {@link ComponentType}. Cannot be null.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @return the registered component.
	 * @throws IllegalArgumentException if pdbccId or type is null,
	 *  or the pdbccId has been registered as a different type,
	 *  or terminal condition is indicated for non-amino-acid component,
	 *  or both N-terminal and C-terminal are true.
	 */
	public static Component register(final String pdbccId, final ComponentType type, 
			final boolean isNTerminal, final boolean isCTerminal) {
		Component comp = of(pdbccId, isNTerminal, isCTerminal);
		if (comp!=null) {
			// already registered
			if (comp.getType()!=type) {
				throw new IllegalArgumentException(pdbccId+" has already been registered" +
						" as a diffent type.");
			}
			return comp;
		}

		comp = new Component(pdbccId, type, isNTerminal, isCTerminal);
		if (isNTerminal) {
			nTerminalAminoAcids.put(pdbccId, comp);
		} else if (isCTerminal) {
			cTerminalAminoAcids.put(pdbccId, comp);
		} else {
			nonTerminalComps.put(pdbccId, comp);
		}
		
		return comp;
	}
	
	/**
	 * Get a non-terminal component.
	 * @param pdbccId Protein Data Bank ID.
	 * @return the Component satisfied.
	 */
	public static Component of(final String pdbccId) {
		return of(pdbccId, false, false);
	}
	
	/**
	 * Get a Component.
	 * @param pdbccId Protein Data Bank ID.
	 * @param isNTerminal true if occurring at N-terminal. false, otherwise.
	 * @param isCTerminal true if occurring at C-terminal. false, otherwise.
	 * @return the Component satisfied.
	 * @throws both N-terminal and C-terminal are specified as true.
	 */
	public static Component of(final String pdbccId, final boolean isNTerminal,
			final boolean isCTerminal) {
		if (isNTerminal&&isCTerminal) {
			throw new IllegalArgumentException("An amino acid can be at" +
			"N-terminal or C-terminal but not both."); //TODO: is this true?
		}
		
		lazyInit();
		
		if (isNTerminal) {
			return nTerminalAminoAcids.get(pdbccId);
		}
		
		if (isCTerminal) {
			return cTerminalAminoAcids.get(pdbccId);
		}
		
		return nonTerminalComps.get(pdbccId);
	}
	
	/**
	 * 
	 * @return informative description.
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getPdbccId());
		sb.append("["+getType().label()+"]");
		if (isCTerminal()) {
			sb.append("(C)");
		} else if (isNTerminal()) {
			sb.append("(N)");
		}
		return sb.toString();
	}
}
