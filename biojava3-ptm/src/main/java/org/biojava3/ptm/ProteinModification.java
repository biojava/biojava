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

package org.biojava3.ptm;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public final class ProteinModification {
	
	/**
	 * Constructor is private, so that we don't
	 * get any stand-alone ProteinModifications. 
	 * ProteinModifications should be obtained from 
	 * getModificationBy... methods. Information about
	 * DataSources can be added with {@link register}.
	 */
	private ProteinModification() {
		throw new AssertionError();
	}
	
	private String id;
	private String pdbccId = null;
	private String pdbccName = null;
	private String residId = null;
	private String residName = null;
	private String psimodId = null;
	private String psimodName = null;
	private String sysName = null;
	
	private boolean isNTerminal = false;
	private boolean isCTerminal = false;
	private String[] components = null;
	private String[] atoms = null;
	
	private String formula;
	private ModificationCategory category;
	private ModificationOccurrenceType occurrenceType;
	
	private static Map<String, ProteinModification> byId = 
		new HashMap<String, ProteinModification>();
	private static Map<String, ProteinModification> byResidId = 
		new HashMap<String, ProteinModification>();
	private static Map<String, ProteinModification> byPsimodId = 
		new HashMap<String, ProteinModification>();
	
	/**
	 * For CHEMICAL_MODIFICATION, one PdbccId should only correspond to 
	 * one instance of ProteinModification. For ATTACHMENT, one PdbccId
	 * could correspond to multiple instances of ProteinModification.
	 */
	private static Map<String, Set<ProteinModification>> byPdbccId = 
		new HashMap<String, Set<ProteinModification>>();
		
	private static Set<ProteinModification> registry = new HashSet<ProteinModification>();

	/**
	 * 
	 * @return modification id.
	 */
	public String id() {
		return id;
	}
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component ID.
	 */
	public String pdbccId() {
		return pdbccId;
	}
	
	/**
	 * 
	 * @return Protein Data Bank Chemical Component name.
	 */
	public String pdbccName() {
		return pdbccName;
	}
	
	/**
	 * 
	 * @return RESID ID.
	 */
	public String residId() {
		return residId;
	}
	
	/**
	 * 
	 * @return RESID name.
	 */
	public String residName() {
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
	public String psimodName() {
		return psimodName;
	}
	
	/**
	 * 
	 * @return Systematic name.
	 */
	public String systematicName() {
		return sysName;
	}
	
	/**
	 * 
	 * @return true if occurring at N-terminal of a protein; false otherwise.
	 */
	public boolean isNTerminal() {
		return isNTerminal;
	}
	
	/**
	 * 
	 * @return true if occurring at C-terminal of a protein; false otherwise.
	 */
	public boolean isCTerminal() {
		return isCTerminal;
	}
	
	/**
	 * Get the components involved.
	 * @return a array of PDBCC ID's of components involved in the modification.
	 *  For CHEMICAL_MODIFICATION, this should only contain the modified amino acid;
	 *  For ATTACHMENT, this should contain two components: the first is the amino acid
	 *  that is attached to, and the second is the attached group.
	 *  For CROSS_OVER, This should contain the two linked amino acids. 
	 */
	public String[] components() {
		return components;
	}
	
	/**
	 * Get the atoms on the components.
	 * @return a array of atoms involved in the modification.
	 *  For CHEMICAL_MODIFICATION, this information is not required and thus could be null;
	 *  For ATTACHMENT, this should contain two atoms: the first is on the amino acid
	 *  that is attached to, and the second is on the attached group.
	 *  For CROSS_OVER, This should contain the atom on the two linked amino acids. 
	 */
	public String[] atoms() {
		return atoms;
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
	public ModificationCategory category() {
		return category;
	}
	
	/**
	 * 
	 * @return the modification occurrence type.
	 */
	public ModificationOccurrenceType occurrenceType() {
		return occurrenceType;
	}
	
	/**
	 * Uses builder pattern to set optional attributes for a ProteinModification. 
	 * For example, this allows you to use the following code:
	 * <pre>
	 * register("0001")
	 *     .residId("AA0406")
	 *     .residName("O-xylosyl-L-serine")
	 *     .components(new String[]{"SER","XYS"})
	 *     .atoms(new String[]{"OG","O1"});
	 * </pre>
	 */
	public static final class Builder {
		private final ProteinModification current;
		
		/**
		 * Create a Builder for a ProteinModification. 
		 * This constructor should only be called by the register method.
		 * @param current the ProteinModification to be modified.
		 */
		private Builder(ProteinModification current) {
			if (current==null)
				throw new IllegalArgumentException("Null argument.");
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
		 * @throws IllegalArgumentException if pdbccId is null or
		 *  it has been set or it is already been registered as 
		 *  CHEMICAL_MODIFICATION or CROSS_OVER.
		 */
		public Builder pdbccId(String pdbccId) {
			if (pdbccId==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.pdbccId!=null) {
				throw new IllegalArgumentException("PDBCC ID has been set.");
			}
			
			current.pdbccId = pdbccId;
			Set<ProteinModification> mods = byPdbccId.get(pdbccId);
			if (mods==null) {
				mods = new HashSet<ProteinModification>();
				byPdbccId.put(pdbccId, mods);
			} else {
				ModificationCategory catByPdbccId = getCategoryByPdbccId(pdbccId);
				if (catByPdbccId!=ModificationCategory.CHEMICAL_MODIFICATION
						|| catByPdbccId!=ModificationCategory.CROSS_OVER) {
					throw new IllegalArgumentException("Only one modification is allowed " +
							"for each PDBCC ID for chemical modifications of cross-over.");
				}
			}
			mods.add(current);
			return this;
		}
		
		/**
		 * Set the Protein Data Bank Chemical Component name.
		 * @param pdbccName Protein Data Bank Chemical Component name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if pdbccName is null or
		 *  it has been set.
		 */
		public Builder pdbccName(String pdbccName) {
			if (pdbccName==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
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
		 * @throws IllegalArgumentException if residId is null or
		 *  it has been set or has been registered by another instance.
		 */
		public Builder residId(String residId) {
			if (residId==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.residId!=null) {
				throw new IllegalArgumentException("RESID ID has been set.");
			}
			
			if (byResidId.containsKey(residId)) {
				// TODO: is this the correct logic?
				throw new IllegalArgumentException(residId+" has been registered.");
			}
			
			current.residId = residId;			
			byResidId.put(residId, current);		
			return this;
		}
		
		/**
		 * Set the RESID name.
		 * @param residName RESID name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if residId is null or
		 *  it has been set.
		 */
		public Builder residName(String residName) {
			if (residName==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
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
		 * @throws IllegalArgumentException if psimodId is null or
		 *  it has been set or has been registered by another instance.
		 */
		public Builder psimodId(String psimodId) {
			if (psimodId==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.psimodId!=null) {
				throw new IllegalArgumentException("PSI-MOD ID has been set.");
			}
			
			if (byResidId.containsKey(psimodId)) {
				// TODO: is this the correct logic?
				throw new IllegalArgumentException(psimodId+" has been registered.");
			}
			
			current.psimodId = psimodId;
			return this;
		}
		
		/**
		 * Set the PSI-MOD name.
		 * @param psimodName PSI-MOD name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if psimodName is null or
		 *  it has been set.
		 */
		public Builder psimodName(String psimodName) {
			if (psimodName==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.psimodName!=null) {
				throw new IllegalArgumentException("PSI-MOD name has been set.");
			}
			
			current.psimodName = psimodName;
			return this;
		}
		
		/**
		 * Set if occurring at N-terminal.
		 * @param isNTerminal true if occurring at N-terminal.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder isNTerminal(boolean isNTerminal) {
			current.isNTerminal = isNTerminal;
			return this;
		}
		
		/**
		 * Set if occurring at C-terminal.
		 * @param isCTerminal true if occurring at C-terminal.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder isCTerminal(boolean isCTerminal) {
			current.isCTerminal = isCTerminal;
			return this;
		}
		
		/**
		 * Set the involved components.
		 * @param components components involved.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if components is null or empty, or
		 *  it has been set, or it has a illegal length. For CHEMICAL_MODIFICATION,
		 *  the length should be 1; otherwise, 2.
		 */
		public Builder components(String[] components) {
			if (components==null || components.length==0) {
				throw new IllegalArgumentException("Null or empty components.");
			}
			
			if (current.components!=null) {
				throw new IllegalArgumentException("Components have been set.");
			}
			
			if (current.category == ModificationCategory.CHEMICAL_MODIFICATION) {
				if (components.length!=1) {
					// TODO: is this necessary?
					throw new IllegalArgumentException("Only one component is allowed for chemical modifications.");
				}
			} else {
				if (components.length!=2) {
					throw new IllegalArgumentException("Two and only two component are required for chemical modifications.");
				}
			}
			
			current.components = components;
			return this;
		}
		
		/**
		 * Set the atoms on the components.
		 * @param atoms atoms on the components.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalStateException if the ProteinModification is a chemical modification.
		 * @throws IllegalArgumentException if atoms is null or length is not 2, or
		 *  it has been set, or the modification is CHEMICAL_MODIFICATION.
		 */
		public Builder atoms(String[] atoms) {
			if (current.category == ModificationCategory.CHEMICAL_MODIFICATION) {
				throw new java.lang.IllegalStateException("Cannot set atoms for chemical modifications.");
			}
			
			if (atoms==null || atoms.length==0) {
				throw new IllegalArgumentException("Null or empty atoms.");
			}
			
			if (current.atoms!=null) {
				throw new IllegalArgumentException("Components have been set.");
			}
			
			if (atoms.length!=2) {
				throw new IllegalArgumentException("Two and only two component are required for chemical modifications.");
			}
			
			current.atoms = atoms;
			return this;
		}
		
		/**
		 * Set the systematic name.
		 * @param sysName systematic name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if sysName is null or
		 *  it has been set.
		 */
		public Builder systematicName(String sysName) {
			if (sysName==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.sysName!=null) {
				throw new IllegalArgumentException("Systematic name has been set.");
			}
			
			current.sysName = sysName;
			return this;
		}
		
		/**
		 * Set the residue formula.
		 * @param formula residue formula.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if formula is null or
		 *  it has been set.
		 */
		public Builder formula(String formula) {
			if (formula==null) {
				throw new IllegalArgumentException("Null argument.");
			}
			
			if (current.formula!=null) {
				throw new IllegalArgumentException("Formula has been set.");
			}
			
			current.formula = formula;
			return this;
		}
	}

	/**
	 * Register a new ProteinModification with (optional) detailed information.
	 * @param id modification id.
	 * @param cat modification category.
	 * @param occType occurrence type.
	 * @return Builder that can be used for adding detailed information.
	 */
	public static Builder register(String id, ModificationCategory cat,
			ModificationOccurrenceType occType) {
		if (id==null || cat==null || occType==null) {
			throw new IllegalArgumentException("Null argument(s)!");
		}
		
		if (byId.containsKey(id)) {
			throw new IllegalArgumentException(id+" has already been registered.");
		}
		
		ProteinModification current = new ProteinModification();
		current.id = id;
		current.category = cat;
		current.occurrenceType = occType;
		
		registry.add(current);
		
		return new Builder(current);
	}

	/**
	 * 
	 * @param residId RESID ID.
	 * @return ProteinModification that has the RESID ID.
	 */
	public static ProteinModification getModificationByResidId(String residId) {
		return byResidId.get(residId);
	}
	/**
	 * 
	 * @param psimodId PSI-MOD ID.
	 * @return ProteinModification that has the PSI-MOD ID.
	 */
	public static ProteinModification getModificationByPsimodId(String psimodId) {
		return byPsimodId.get(psimodId);
	}
	
	/**
	 * 
	 * @param pdbccId
	 * @return
	 */
	public static ProteinModification getChemicalModificationByPdbccId(String pdbccId) {
		ModificationCategory cat = getCategoryByPdbccId(pdbccId);		
		if (cat==null) {
			return null;
		}
		
		if (cat!=ModificationCategory.CHEMICAL_MODIFICATION) {
			// TODO: throw a exception or just return null?
			throw new IllegalArgumentException(pdbccId+" is not a chemical modification.");
		}
		
		return byPdbccId.get(byPdbccId).iterator().next();
	}
	
	public static ModificationCategory getCategoryByPdbccId(String pdbccId) {
		Set<ProteinModification> mods = byPdbccId.get(byPdbccId);
		if (mods==null||mods.isEmpty()) {
			return null;
		}
		
		if (mods.size()>1) {
			return ModificationCategory.ATTACHMENT;
		}
		
		return mods.iterator().next().category();
	}
}
