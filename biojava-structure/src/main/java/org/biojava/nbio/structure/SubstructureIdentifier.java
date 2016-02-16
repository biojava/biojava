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
 * Created on December 19, 2013
 * Author: Douglas Myers-Turnbull
 */

package org.biojava.nbio.structure;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is the canonical way to identify a part of a structure.
 * 
 * <p>The current syntax allows the specification of a set of residues from
 * the first model of a structure. Future versions may be extended to represent
 * additional properties.
 * 
 * <p>Identifiers should adhere to the following specification, although some
 * additional forms may be tolerated where unambiguous for backwards compatibility.
 * <pre>
 * 		name          := pdbID
 * 		               | pdbID '.' chainID
 * 		               | pdbID '.' range
 * 		range         := range (',' range)?
 * 		               | chainID
 * 		               | chainID '_' resNum '-' resNum
 * 		pdbID         := [0-9][a-zA-Z0-9]{3}
 * 		chainID       := [a-zA-Z0-9]+
 * 		resNum        := [-+]?[0-9]+[A-Za-z]?
 * </pre>
 * For example:
 * <pre>
 * 		1TIM                            #whole structure
 * 		1tim                            #same as above
 * 		4HHB.C                          #single chain
 * 		3AA0.A,B                        #two chains
 * 		4GCR.A_1-40                     #substructure
 *      3iek.A_17-28,A_56-294,A_320-377 #substructure of 3 disjoint parts
 * </pre>
 * More options may be added to the specification at a future time.

 * @author dmyersturnbull
 * @author Spencer Bliven
 */
public class SubstructureIdentifier implements Serializable, StructureIdentifier {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(SubstructureIdentifier.class);

	private final String pdbId;
	private final List<ResidueRange> ranges;
	
	/**
	 * Create a new identifier from a string.
	 * @param id
	 */
	public SubstructureIdentifier(String id) {
		String[] idRange = id.split("\\.");
		if(1 > idRange.length || idRange.length > 2 ) {
			throw new IllegalArgumentException(String.format("Malformed %s: %s",getClass().getSimpleName(),id));
		}
		this.pdbId = idRange[0];
		if(this.pdbId.length() != 4) {
			// Changed from Exception to a warning to support files and stuff -sbliven 2015/01/22
			logger.warn(String.format("Unrecognized PDB code %s",this.pdbId));
		}
		
		if( idRange.length == 2) {
			String rangeStr = idRange[1].trim();

			this.ranges = ResidueRange.parseMultiple(rangeStr);
		} else {
			this.ranges = new LinkedList<ResidueRange>();
		}
	}

	/**
	 * Create a new identifier based on a set of ranges.
	 * 
	 * If ranges is empty, includes all residues.
	 * @param pdbId
	 * @param ranges
	 */
	public SubstructureIdentifier(String pdbId, List<ResidueRange> ranges) {
		if(ranges == null) {
			throw new NullPointerException("Null ranges list");
		}
		this.pdbId = pdbId;
		this.ranges = ranges;
	}

	@Override
	public String toString() {
		return getIdentifier();
	}

	/**
	 * Get the String form of this identifier.
	 * 
	 * This provides the canonical form for a StructureIdentifier and has
	 * all the information needed to recreate a particular substructure.
	 *
	 * Example: 3iek.A_17-28,A_56-294
	 * @return The String form of this identifier
	 */
	@Override
	public String getIdentifier() {
		if (ranges.isEmpty()) return pdbId;
		return pdbId + "." + ResidueRange.toString(ranges);
	}

	public String getPdbId() {
		return pdbId;
	}

	public List<ResidueRange> getResidueRanges() {
		return ranges;
	}

	/**
	 * Return itself. SubstructureIdentifiers are canonical!
	 */
	@Override
	public SubstructureIdentifier toCanonical() {
		return this;
	}

	/**
	 * Takes a complete structure as input and reduces it to residues present in
	 * the specified ranges
	 * 
	 * <p>The returned structure will be a shallow copy of the input, with shared
	 * Chains, Residues, etc.
	 * @param input A full structure, e.g. as loaded from the PDB. The structure
	 * ID should match that returned by getPdbId().
	 * @return 
	 * @throws StructureException 
	 * @see StructureTools#getReducedStructure(Structure, String)
	 */
	@Override
	public Structure reduce(Structure s) throws StructureException {
		// Follows StructureImpl.clone()

		// Create new structure & copy basic properties
		Structure newS = new StructureImpl();

		newS.setPDBCode(s.getPDBCode());
		newS.setPDBHeader(s.getPDBHeader());
		newS.setName(this.toString());
		newS.setDBRefs(s.getDBRefs());
		newS.setBiologicalAssembly(s.isBiologicalAssembly());
		newS.getPDBHeader().setDescription(
				"sub-range " + ranges + " of "  + newS.getPDBCode() + " "
						+ s.getPDBHeader().getDescription());
		newS.setCrystallographicInfo(s.getCrystallographicInfo());
		// TODO The following should be only copied for atoms which are present in the range.
		newS.setCompounds(s.getCompounds());
		newS.setConnections(s.getConnections());
		newS.setSSBonds(s.getSSBonds());
		newS.setSites(s.getSites());

		newS.setStructureIdentifier(this);

		for( int modelNr=0;modelNr<s.nrModels();modelNr++) {
			String prevChainId = null;
			
			
			// Construct new model
			newS.addModel(new ArrayList<Chain>());
			
			if(getResidueRanges().isEmpty()) {
				// Include all residues
				newS.setCompounds(s.getCompounds());
				newS.setConnections(s.getConnections());
				newS.setSSBonds(s.getSSBonds());
				newS.setSites(s.getSites());

				newS.setModel(modelNr, s.getModel(modelNr));
			} else {
				// Restrict residues
				for( ResidueRange range: getResidueRanges()) {

					String chainId = range.getChainId();
					ResidueNumber pdbresnum1 = range.getStart();
					ResidueNumber pdbresnum2 = range.getEnd();

					Chain chain;
					if(chainId.equals("_") ) {
						// Handle special case of "_" chain for single-chain proteins
						chain = s.getChain(modelNr,0);
						if(pdbresnum1 != null)
							pdbresnum1.setChainId(chain.getChainID());
						if(pdbresnum2 != null)
							pdbresnum2.setChainId(chain.getChainID());

						if(s.size() != 1) {
							// SCOP 1.71 uses this for some proteins with multiple chains
							// Print a warning in this ambiguous case
							logger.warn("Multiple possible chains match '_'. Using chain {}",chain.getChainID());
						}
					} else {
						// Explicit chain
						try {
							chain = s.getChainByPDB(chainId,modelNr);
						} catch(StructureException e) {
							// Chain not found
							// Maybe it was a chain index, masquerading as a chainId?
							try {
								int chainNum = Integer.parseInt(chainId);
								try {
									chain = s.getChain(modelNr, chainNum);
									logger.warn("No chain found for {}. Interpretting it as an index, using chain {} instead",chainId,chain.getChainID());
								} catch(Exception e2) { //we don't care what gets thrown here -sbliven
									throw e; // Nope, not an index. Throw the original exception
								}
							} catch(NumberFormatException e3) {
								// Not an index. Throw the original exception
								throw e;
							}
						}
					}

					List<Group> groups;
					if(pdbresnum1 == null && pdbresnum2 == null) {
						groups = chain.getAtomGroups();
					} else {
//						// Trim extra residues off the range
//						Atom[] allAtoms = StructureTools.getRepresentativeAtomArray(chain);
//						AtomPositionMap map = new AtomPositionMap(allAtoms);
//						ResidueRange trimmed = map.trimToValidResidues(
//								new ResidueRange(chain.getChainID(),
//										pdbresnum1, pdbresnum2));
//						if (trimmed != null) {
//							pdbresnum1 = trimmed.getStart();
//							pdbresnum2 = trimmed.getEnd();
//						}
						groups = Arrays.asList(chain.getGroupsByPDB(pdbresnum1, pdbresnum2));
					}

					// Create new chain, if needed
					Chain c = null;
					if ( prevChainId == null) {
						// first chain...
						c = new ChainImpl();
						c.setChainID(chain.getChainID());
						newS.addChain(c,modelNr);
					} else if ( prevChainId.equals(chain.getChainID())) {
						c = newS.getChainByPDB(prevChainId,modelNr);

					} else {
						try {
							c = newS.getChainByPDB(chain.getChainID(),modelNr);
						} catch (StructureException e){
							// chain not in structure yet...
							c = new ChainImpl();
							c.setChainID(chain.getChainID());
							newS.addChain(c,modelNr);
						}
					}

					// add the groups to the chain:
					for ( Group g: groups) {
						c.addGroup(g);
					}

					prevChainId = c.getChainID();
				} // end range
			}
		} // end modelNr

		return newS;
	}

	/**
	 * Loads the complete structure based on {@link #getPdbId()}.
	 * 
	 * @param AtomCache A source of structures
	 * @return A Structure containing at least the atoms identified by this,
	 *  or null if no PDB ID is set
	 * @throws StructureException For errors loading and parsing the structure
	 * @throws IOException Errors reading the structure from disk
	 */
	@Override
	public Structure loadStructure(AtomCache cache) throws IOException, StructureException {
		String pdb = getPdbId();
		if(pdb == null)
			return null;
		return cache.getStructureForPdbId(pdb);
	}

}
