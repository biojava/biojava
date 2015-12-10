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
 */
package org.biojava.nbio.structure.align.client;


import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.PassthroughIdentifier;
import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.domain.PDPDomain;
import org.biojava.nbio.structure.domain.PDPProvider;
import org.biojava.nbio.structure.domain.RemotePDPProvider;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.biojava.nbio.structure.scop.ScopFactory;


/** A utility class that makes working with names of structures, domains and ranges easier.
 * 
 * @see #getName the name. e.g. 4hhb, 4hhb.A, d4hhba_, PDP:4HHBAa etc.
 */
public class StructureName implements Comparable<StructureName>, Serializable, StructureIdentifier {
	
	//private static final Logger logger = LoggerFactory.getLogger(StructureName.class);

	private static final long serialVersionUID = 4021229518711762957L;
	protected String name;
	protected String pdbId;
	protected String chainId;

	private static final Pattern cathPattern = Pattern.compile("^([0-9][a-z0-9]{3})(\\w)([0-9]{2})$",Pattern.CASE_INSENSITIVE);
	// More specific than AtomCache.scopIDregex
	private static final Pattern scopPattern = Pattern.compile("^d([0-9][a-z0-9]{3})(\\w|\\.)(\\w)$",Pattern.CASE_INSENSITIVE);

	private enum Source {
		PDB,
		SCOP,
		PDP,
		CATH,
		URL,
		FILE,
	};

	private Source mySource = null; 

	// cache for realize() method
	private StructureIdentifier realized = null;
	
	/**
	 * Create a new StructureName from the given identifier, which may be a 
	 * domain name, a substructure identifier, etc.
	 * 
	 * The source and PDB-Id are extracted at compile time, but fully
	 * interpreting the ID, which may require additional parsing or remote
	 * calls, is done lazily.
	 * @param name An identifier string
	 * @throws IllegalArgumentException if the name has a recognizable source but is semantically invalid
	 */
	public StructureName(String name){
		this.name = name;

		init();//sets pdbId and mySource
	}


	/**
	 * Tries to determine the source and pdbId without fully realizing the identifier,
	 * which could require I/O depending on the source
	 * @throws IllegalArgumentException if the source is recognizable but invalid
	 */
	private void init(){
		// SCOP domain
		Matcher matcher = scopPattern.matcher(name);
		if ( matcher.matches() ) {
			mySource = Source.SCOP;
			pdbId = matcher.group(1);
			chainId = matcher.group(2);
			return;
		}
		// PDP
		if ( name.startsWith(AtomCache.PDP_DOMAIN_IDENTIFIER)){
			// starts with PDP:
			// eg: PDP:3LGFAa
			mySource = Source.PDP;
			matcher = PDPDomain.PDP_NAME_PATTERN.matcher(name);
			if(! matcher.matches() ) {
				throw new IllegalArgumentException("Malformed PDP domain name");
			}
			pdbId = matcher.group(1);
			chainId = matcher.group(2);
			return;
		}
		// CATH
		matcher = cathPattern.matcher(name);
		if ( matcher.matches() ){
			mySource = Source.CATH;
			pdbId = matcher.group(1);
			chainId = matcher.group(2);
			return;
		}
		// URL
		try {
			new URL(name);
			mySource = Source.URL;
			pdbId = null;
			chainId = null;
			return;
		} catch(MalformedURLException e) {}
		// File
		// TODO detect non-existent files
		File file = new File(FileDownloadUtils.expandUserHome(name));
		if( file.exists() ) {
			mySource = Source.FILE;
			pdbId = null;
			chainId = null;
			return;
		}

		// Default to PDB
		mySource = Source.PDB;
		SubstructureIdentifier si = new SubstructureIdentifier(getIdentifier());
		realized = si; // Safe to realize immediately

		pdbId = si.getPdbId();
		// Set chainId if unique
		Set<String> chains = getChainIds(si);
		if(chains.size() == 1) {
			this.chainId = chains.iterator().next();
		} else if(chains.size() > 1) {
			this.chainId = ".";
		} else {
			this.chainId = null;
		}
	}

	private static Set<String> getChainIds(SubstructureIdentifier si) {
		Set<String> chains = new TreeSet<String>();
		List<ResidueRange> ranges = si.getResidueRanges();
		for(ResidueRange range : ranges) {
			String chain = range.getChainId();
			if(chain != null) {
				chains.add(chain);
			}
		}
		return chains;
	}

	/**
	 * Get the PDB ID for this name, if any.
	 * 
	 * Equivalent to {@link SubstructureIdentifier#getPdbId()
	 * toCanonical().getPdbId()}
	 * @return The upper-case PDB Name, or null if not applicable
	 * @throws StructureException Wraps errors which occur when converting to canonical form
	 */
	public String getPdbId() throws StructureException {
		if( pdbId == null) {
			pdbId = toCanonical().getPdbId();
		}
		return pdbId.toUpperCase();
	}

	/**
	 * Gets the chain ID, for structures where it is unique and well-defined.
	 * May return '.' for multi-chain ranges, '_' for wildcard chains, or
	 * null if the information is unavailable.
	 * 
	 * <p>This method should only be used casually. For precise chainIds, it
	 * is better to use {@link #toCanonical()} and iterate through the
	 * residue ranges.
	 * @return
	 */
	public String getChainId() {
		return chainId;
	}
	/**
	 * 
	 * @return the identifier string
	 * @deprecated use {@link #getIdentifier()}
	 */
	@Deprecated
	public String getName(){

		return getIdentifier();
	}

	@Override
	public String getIdentifier() {
		return name;
	}

	@Override
	public String toString(){

		StringBuilder s = new StringBuilder();

		s.append(name);
		s.append("[");

		try {
			String id = getPdbId();
			s.append("PDB ID: ");
			s.append(id);
			s.append(",");
		} catch (StructureException e) {
			// ignore errors
		}
		s.append("Source: ");
		s.append(mySource);
		s.append("]");

		return s.toString();

	}


	public boolean isScopName() {
		return mySource == Source.SCOP;
	}

	public boolean isPDPDomain(){
		return mySource == Source.PDP;
	}

	public boolean isCathID(){
		return mySource == Source.CATH;
	}

	public boolean isPdbId(){
		return mySource == Source.PDB;
	}

	public boolean isURL() {
		return mySource == Source.URL;
	}

	public boolean isFile() {
		return mySource == Source.FILE;
	}

	/**
	 * 
	 * @return
	 * @throws StructureException Wraps exceptions that may be thrown by
	 *  individual implementations. For example, a SCOP identifier may require
	 *  that the domain definitions be available for download.
	 */
	private StructureIdentifier realize() throws StructureException {
		if( realized == null ) {

			switch(mySource) {
			case CATH:
				realized = CathFactory.getCathDatabase().getDescriptionByCathId(getIdentifier());
				break;
			case SCOP:
				realized = ScopFactory.getSCOP().getDomainByScopID(getIdentifier());
				break;
			case FILE:
			case URL:
				realized = new PassthroughIdentifier(name);
				break;
			case PDP:
				//TODO -sbliven 2015-01-28
				try {
					PDPProvider provider = new RemotePDPProvider(false);
					realized = provider.getPDPDomain(name);
				} catch (IOException e) {
					// This is really bad, but the SCOP and CATH factories do it internally too -sbliven
					throw new StructureException("Unable to fetch PDP domain "+name, e);
				}
				break;
			case PDB:
			default:
				realized = new SubstructureIdentifier(getIdentifier());
				break;
			}
		}
		return realized;
	}

	@Override
	public SubstructureIdentifier toCanonical() throws StructureException {
		return realize().toCanonical();
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return realize().reduce(input);
	}
	
	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return realize().loadStructure(cache);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		StructureName other = (StructureName) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	/**
	 * Orders identifiers lexicographically by PDB ID and then full Identifier
	 */
	@Override
	public int compareTo(StructureName o) {
		if ( this.equals(o))
			return 0;

		String pdb1 = null;
		String pdb2 = null;
		try {
			pdb1 = this.getPdbId();
		} catch (StructureException e) {}
		try {
			pdb2 = this.getPdbId();
		} catch (StructureException e) {}

		int comp = 0;

		// Sort those with PDBIDs before those without
		if( pdb1 == null ) {
			if( pdb2 != null) {
				return 1; // this > o
			}
			// both null
		} else if( pdb2 == null){
			return -1; // this < o
		} else {
			// neither null
			comp = pdb1.compareTo(pdb2);
		}
		if( comp != 0 ) {
			return comp;
		}

		// break tie with full identifiers
		pdb1 = this.getIdentifier();
		pdb2 = o.getIdentifier();

		// Throws NPE for nulls
		return pdb1.compareTo(pdb2);
	}
}
