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

/**
 * contains information about a certain ProteinModification.
 * The ProteinModification class uses the extensible enum pattern.
 * You can't instantiate ProteinModifications directly, instead 
 * you have to use one of the getBy... methods.
 */
public final class ProteinModification {
	
	/**
	 * Constructor is private, so that we don't
	 * get any stand-alone ProteinModifications. 
	 * ProteinModifications should be obtained from 
	 * getBy... methods. Information about
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
	private String formula = null;
	private String description = null;
	
	private boolean isNTerminal = false;
	private boolean isCTerminal = false;
	private String[] components = null;
	private String[][] atoms = null;
	
	private ModificationCategory category;
	private ModificationOccurrenceType occurrenceType;
	
	private final static Map<String, ProteinModification> byId = 
		new HashMap<String, ProteinModification>();
	private final static Map<String, ProteinModification> byResidId = 
		new HashMap<String, ProteinModification>();
	private final static Map<String, ProteinModification> byPsimodId = 
		new HashMap<String, ProteinModification>();
	private final static Map<String, Set<ProteinModification>> byPdbccId = 
		new HashMap<String, Set<ProteinModification>>();		
	private static Set<ProteinModification> registry = 
		new HashSet<ProteinModification>();

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
	 * @return Description.
	 */
	public String description() {
		return description;
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
	 *  For CROSS_OVER, This should contain the linked amino acids and other 
	 *  chemical components. 
	 */
	public String[] components() {
		return components;
	}
	
	/**
	 * Get the atoms on the components.
	 * @return a matrix of atoms involved in the modification. A non-null element in the
	 *  matrix indicates a link between the corresponding two elements.
	 *  For CHEMICAL_MODIFICATION, this information is not required and thus could be null;
	 *  For ATTACHMENT, this should be a 2x2 matrix: element (1,2) represents the atom on 
	 *  the amino acid that is attached to, and element (2,1) represents the atom on the 
	 *  attached group.
	 *  For CROSS_OVER, this should be a c by c matrix, where c is the size of components.
	 */
	public String[][] atoms() {
		return atoms;
	}
	
	/**
	 * 
	 * @return formula of the modified residue.
	 */
	public String formula() {
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
	 * ProteinModification.
	 *     .register("0001", Modification.ATTACHMENT, ModificationOccurrenceType.NATURAL)
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
		private Builder(final ProteinModification current) {
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
		 * @throws IllegalArgumentException if pdbccId has been set.
		 */
		public Builder pdbccId(final String pdbccId) {			
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
		public Builder pdbccName(final String pdbccName) {			
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
		 * @throws IllegalArgumentException if residIdhas been set 
		 * or has been registered by another instance.
		 */
		public Builder residId(final String residId) {
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
		 * @throws IllegalArgumentException if residId has been set.
		 */
		public Builder residName(final String residName) {
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
		 * @throws IllegalArgumentException if psimodId has been set
		 *  or has been registered by another instance.
		 */
		public Builder psimodId(final String psimodId) {
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
		 * @throws IllegalArgumentException if psimodName has been set.
		 */
		public Builder psimodName(final String psimodName) {
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
		public Builder isNTerminal(final boolean isNTerminal) {
			current.isNTerminal = isNTerminal;
			return this;
		}
		
		/**
		 * Set if occurring at C-terminal.
		 * @param isCTerminal true if occurring at C-terminal.
		 * @return the same Builder object so you can chain setters.
		 */
		public Builder isCTerminal(final boolean isCTerminal) {
			current.isCTerminal = isCTerminal;
			return this;
		}
		
		/**
		 * Set one component. The method is for ModifiedResidue only.
		 * @param component component.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if component is null.
		 */
		public Builder component(final String component) {
			if (component==null) {
				throw new IllegalArgumentException("Null component.");
			}
			
			return componentsAndAtoms(new String[]{component}, null);
		}
		
		/**
		 * Set components and atoms for modification involving two components,
		 * e.g. attachments.
		 * @param component1 the first involved component (amino acid for attachment).
		 * @param atom1 atom on the first component.
		 * @param component2 the second involved component (attached group for attachment).
		 * @param atom2 atom on the second component.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if any of the arguments is null.
		 */
		public Builder componentsAndAtoms(final String component1,
				final String atom1, final String component2, final String atom2) {
			if (component1==null || atom1==null || component2==null || atom2==null) {
				throw new IllegalArgumentException("Null argument(s).");
			}
			
			return componentsAndAtoms(new String[]{component1, component2},
					new String[][]{new String[]{null,atom1},new String[]{atom2,null}});
		}
		
		/**
		 * Set the involved components and atoms.
		 * @param components components involved.
		 * @param atoms a matrix of atoms involved in the modification. A non-null element 
		 *  in the matrix indicates a link between the corresponding two elements.
		 *  For CHEMICAL_MODIFICATION, this information is not required and thus could be null;
		 *  For ATTACHMENT, this should be a 2x2 matrix: element (1,2) represents the atom on 
		 *  the amino acid that is attached to, and element (2,1) represents the atom on the 
		 *  attached group.
		 *  For CROSS_OVER, this should be a c by c matrix, where c is the size of components.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if components or atoms has been set,
		 *  or atom matrix dimension is improper.
		 */
		public Builder componentsAndAtoms(final String[] components,
				final String[][] atoms) {
			if (current.components!=null) {
				throw new IllegalArgumentException("Components have been set.");
			}
			
			if (current.atoms!=null) {
				throw new IllegalArgumentException("atoms have been set.");
			}
			
			if (components==null || components.length==0) {
				return this;
			}			
			
			if (atoms!=null) {
				if (components.length!=atoms.length || components.length!=atoms[0].length) {
					throw new IllegalArgumentException("The matrix of atoms must has the same" +
							"number of elements as components in both dimensions.");
				}
			}
			
			current.components = components;
			current.atoms = atoms;
			
			return this;
		}
		
		/**
		 * Set the systematic name.
		 * @param sysName systematic name.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if sysName has been set.
		 */
		public Builder systematicName(final String sysName) {			
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
		public Builder description(final String description) {
			if (current.description!=null) {
				throw new IllegalArgumentException("Description has been set.");
			}
			
			current.description = description;
			
			return this;
		}
		
		/**
		 * Set the residue formula.
		 * @param formula residue formula.
		 * @return the same Builder object so you can chain setters.
		 * @throws IllegalArgumentException if formula has been set.
		 */
		public Builder formula(final String formula) {
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
	public static Builder register(final String id, final ModificationCategory cat,
			final ModificationOccurrenceType occType) {
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
	 * @param id modification ID.
	 * @return ProteinModification that has the corresponding ID.
	 */
	public static ProteinModification getById(final String id) {
		return byId.get(id);
	}

	/**
	 * 
	 * @param residId RESID ID.
	 * @return ProteinModification that has the RESID ID.
	 */
	public static ProteinModification getByResidId(final String residId) {
		return byResidId.get(residId);
	}
	/**
	 * 
	 * @param psimodId PSI-MOD ID.
	 * @return ProteinModification that has the PSI-MOD ID.
	 */
	public static ProteinModification getByPsimodId(final String psimodId) {
		return byPsimodId.get(psimodId);
	}
	
	/**
	 * 
	 * @param pdbccId Protein Data Bank Chemical Component ID.
	 * @return chemical modifications that have the PDBCC ID.
	 */
	public static Set<ProteinModification> getByPdbccId(final String pdbccId) {
		return byPdbccId.get(byPdbccId);
	}
}
