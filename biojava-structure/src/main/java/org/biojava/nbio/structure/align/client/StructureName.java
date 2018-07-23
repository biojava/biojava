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
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.BioAssemblyIdentifier;
import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.URLIdentifier;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.domain.PDPDomain;
import org.biojava.nbio.structure.domain.PDPProvider;
import org.biojava.nbio.structure.domain.RemotePDPProvider;
import org.biojava.nbio.structure.ecod.EcodFactory;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * A utility class that makes working with names of structures, domains and ranges easier.
 *
 * Accepts a wide range of identifier formats, including {@link ScopDomain},
 * {@link CathDomain}, PDP domains, and {@link SubstructureIdentifier} residue
 * ranges.
 *
 * Where possible, data is extracted from the input string. Otherwise, range
 * information may be loaded from one of the factory classes:
 * {@link CathFactory},{@link ScopFactory}, etc.
 *
 * @see #getName the name. e.g. 4hhb, 4hhb.A, d4hhba_, PDP:4HHBAa etc.
 */

public class StructureName implements Comparable<StructureName>, Serializable, StructureIdentifier {
	private static final long serialVersionUID = 4021229518711762957L;
	private static final Logger logger = LoggerFactory.getLogger(StructureName.class);

	protected String name;
	protected String pdbId;
	protected String chainName;

	private static final Pattern cathPattern = Pattern.compile("^(?:CATH:)?([0-9][a-z0-9]{3})(\\w)([0-9]{2})$",Pattern.CASE_INSENSITIVE);
	// ds046__ is a special case with no PDB entry
	private static final Pattern scopPattern = Pattern.compile("^(?:SCOP:)?d([0-9][a-z0-9]{3}|s046)(\\w|\\.)(\\w)$",Pattern.CASE_INSENSITIVE);
	// ECOD chains and domains can't be automatically distinguished. Ex: e3j9zS13 is chain 'S1', e1wz2B14 is chain 'B'
	private static final Pattern ecodPattern = Pattern.compile("^(?:ECOD:)?e([0-9][a-z0-9]{3})(?:\\w|\\.)\\w+$",Pattern.CASE_INSENSITIVE);

	// Names are automatically used as prefixes
	public enum Source {
		PDB,
		SCOP,
		PDP,
		CATH,
		URL,
		FILE,
		ECOD,
		BIO,
	};

	private Source mySource = null;

	// cache for getBaseIdentifier() method
	private StructureIdentifier base = null;

	/**
	 * Create a new StructureName from the given identifier, which may be a
	 * domain name, a substructure identifier, etc.
	 * <p>
	 * The source and PDB-Id are extracted at compile time, but fully
	 * interpreting the ID, which may require additional parsing or remote
	 * calls, is done lazily.
	 * <p>
	 * The following sources are supported. Any may be prefixed by the source
	 * name followed by a colon (e.g. PDB:4HHB). In this case, that source will be used
	 * unequivocally. If no source is specified, StructureName will make a
	 * (usually reliable) guess as to which source was intended.
	 * <ul>
	 * <li><b>PDB</b>PDB identifier, optionally followed by chain and/or residue
	 *     ranges. Internally represented by a {@link SubstructureIdentifier};
	 *     see that class for the full format specification.
	 *     Examples: 4hhb, 4hhb.A, 4hhb.A:1-50.
	 * <li><b>SCOP</b> SCOP domain (or SCOPe, depending on the
	 *     {@link ScopFactory#getSCOP()} version). Example: d1h6w.2
	 * <li><b>PDP</b> Protein Domain Parser domain. PDP domains are not guessed,
	 *     making the PDP: prefix obligatory. Example: PDP:4HHBAa
	 * <li><b>CATH</b> Cath domains. Example: 1qvrC03
	 * <li><b>URL</b> Arbitrary URLs. Most common protocols are handled,
	 *     including http://, ftp://, and file://. Some parsing information can
	 *     be passed as custom query parameters. Example:
	 *     http://www.rcsb.org/pdb/files/1B8G.pdb.gz
	 * <li><b>FILE</b> A file path. Supports relative paths and expands ~ to
	 *     the user's home directory. Only existing files will be automatically
	 *     detected; to refer to a potentially not-yet existing file, prepend
	 *     the prefix. Internally represented as a {@link URLIdentifier}
	 *     after path expansion. Example: ~/custom_protein.pdb
	 * <li><b>ECOD</b> ECOD domain. Example: e1lyw.1
	 * <li><b>BIO</b> Biological assembly. These are not guessed, making
	 *     the BIO: prefix obligatory. Example: BIO:2ehz:1
	 * </ul>
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

		// First try identifying a prefix
		String[] prefix = name.split(":", 2);
		mySource = null;
		if(prefix.length > 1) {
			// Match Source prefixes
			String suffix = prefix[1];
			try {
				mySource = Source.valueOf(prefix[0].toUpperCase());
			} catch( IllegalArgumentException e ) {
				// unrecognized prefix; fall back on guessing
				mySource = null;
			}
			if(mySource != null) {
				switch( mySource) {
				case SCOP:
					if( ! initFromScop(suffix) )
						throw new IllegalArgumentException("Malformed SCOP domain name:"+suffix);
					return;
				case PDP:
					if( ! initFromPDP(name) )
						throw new IllegalArgumentException("Malformed PDP domain name:"+suffix);
					return;
				case CATH:
					if( ! initFromCATH(suffix) )
						throw new IllegalArgumentException("Malformed CATH domain name:"+suffix);
					return;
				case BIO:
					if( ! initFromBIO(name) )
						throw new IllegalArgumentException("Malformed BIO name:"+suffix);
					return;
				case ECOD:
					if( ! initFromECOD(suffix) )
						throw new IllegalArgumentException("Malformed ECOD domain name:"+suffix);
					return;
				case PDB:
					if( ! initFromPDB(suffix) )
						throw new IllegalArgumentException("Malformed PDB specification:"+suffix);
					return;
				case FILE:
					// Treat file:/ prefixes as URLs
					if( ! suffix.startsWith("/")) {
						// Otherwise, treat as file
						initFromFile();
						return;
					}
					// fall through to URL case
				case URL:
					if( ! initFromURL(name))
						throw new IllegalArgumentException("Malformed URL specification:"+suffix);
					return;
				default:
					throw new IllegalStateException("Unimplemented Source "+mySource);
				}
			}
		}

		// No known prefix, so revert to guessing

		// First guess regex-based identifiers
		// SCOP domain
		if( initFromScop(name) )
			return;
		// CATH
		if( initFromCATH(name) )
			return;
		// ECOD
		if( initFromECOD(name) )
			return;
		// Never guess BIO or PDP

		// URL
		if( initFromURL(name) )
			return;

		// Guess FILE based on file existence
		File file = new File(FileDownloadUtils.expandUserHome(name));
		if( file.canRead() && !file.isDirectory() ) {
			// an attempt to mitigate issue #398. It doesn't fix it but it catches the most common case of passing a pdb id and finding a file in working dir matching it
			if (name.matches("\\d\\w\\w\\w")) {
				// the plain pdb id case, this is unlikely to be what the user wants: let's let it through but warn about it
				logger.warn("Provided 4-letter structure name '{}' matches "
						+ "file name in directory {}. Will read structure "
						+ "data from file {} and not consider the name as a "
						+ "structure identifier. If this is not what you "
						+ "want, use 'FILE:{}'",
						name, file.getAbsoluteFile().getParent(),
						file.getAbsolutePath(), name);
			} else {
				logger.info("Provided structure name '{}' matches "
						+ "file name in directory {}. Will read structure "
						+ "data from file {}.",
						name, file.getAbsoluteFile().getParent(),
						file.getAbsolutePath());
			}

			initFromFile();
			return;
		}

		// Default to PDB
		initFromPDB( name );
	}

	private boolean initFromScop(String name) {
		Matcher matcher = scopPattern.matcher(name);
		if ( matcher.matches() ) {
			mySource = Source.SCOP;
			pdbId = matcher.group(1).toUpperCase();
			chainName = matcher.group(2);
			return true;
		}
		return false;
	}
	private boolean initFromPDP(String name) {
		Matcher matcher = PDPDomain.PDP_NAME_PATTERN.matcher(name);
		if( matcher.matches() ) {
			pdbId = matcher.group(1).toUpperCase();
			chainName = matcher.group(2);
			return true;
		}
		return false;
	}
	private boolean initFromCATH(String name) {
		Matcher matcher = cathPattern.matcher(name);
		if ( matcher.matches() ){
			mySource = Source.CATH;
			pdbId = matcher.group(1).toUpperCase();
			chainName = matcher.group(2);
			return true;
		}
		return false;
	}
	private boolean initFromECOD(String name) {
		Matcher matcher = ecodPattern.matcher(name);
		if ( matcher.matches() ){
			mySource = Source.ECOD;
			pdbId = matcher.group(1).toUpperCase();
			chainName = null;
			return true;
		}
		return false;
	}
	private boolean initFromBIO(String name) {
		Matcher matcher = BioAssemblyIdentifier.BIO_NAME_PATTERN.matcher(name);
		if( matcher.matches() ) {
			pdbId = matcher.group(1).toUpperCase();
			return true;
		}
		return false;
	}
	private boolean initFromPDB(String suffix) {
		mySource = Source.PDB;
		SubstructureIdentifier si = new SubstructureIdentifier(suffix);
		base = si; // Safe to realize immediately

		pdbId = si.getPdbId();
		// Set chainName if unique
		Set<String> chains = getChainNames(si);
		if(chains.size() == 1) {
			this.chainName = chains.iterator().next();
		} else if(chains.size() > 1) {
			this.chainName = ".";
		} else {
			this.chainName = null;
		}
		return true;
	}
	private boolean initFromURL(String suffix) {
		try {
			URL url = new URL(suffix);
			String path = url.getPath();
			mySource = Source.URL;
			pdbId = URLIdentifier.guessPDBID( path.substring(path.lastIndexOf('/')+1) );
			chainName = null; // Don't bother checking query params here
			return true;
		} catch(MalformedURLException e) {
			return false;
		}
	}
	private boolean initFromFile() {
		mySource = Source.FILE;
		pdbId = null;
		chainName = null;
		return true;
	}

	private static Set<String> getChainNames(SubstructureIdentifier si) {
		Set<String> chains = new TreeSet<String>();
		List<ResidueRange> ranges = si.getResidueRanges();
		for(ResidueRange range : ranges) {
			String chainName = range.getChainName();
			if(chainName != null) {
				chains.add(chainName);
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
		return pdbId;
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
		return chainName;
	}
	
	/**
	 * Get the original form of the identifier
	 */
	@Override
	public String getIdentifier() {
		return name;
	}

	@Override
	public String toString(){

		return name;
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

	/**
	 * Indicates that the identifier was determined to correspond to a file.
	 * Note that some file identifiers may also be valid URLs; in that case,
	 * the URL source is preferred.
	 * @return
	 */
	public boolean isFile() {
		return mySource == Source.FILE;
	}

	public boolean isEcodDomain() {
		return mySource == Source.ECOD;
	}

	public boolean isBioAssembly() {
		return mySource == Source.BIO;
	}

	public Source getSource() {
		return mySource;
	}

	/**
	 * StructureName wraps another StructureIdentifier. The type of the base
	 * identifier depends on the {@link #getSource() source}. Most StructureName
	 * methods deligate to the base identifier.
	 *
	 * <p>It is possible that future versions of StructureName might change the
	 * return type. Except for some specialized uses, it is probably better
	 * to create the correct type of identifier directly, rather than creating
	 * a StructureName and casting the result of this method.
	 * @return A Str
	 * @throws StructureException Wraps exceptions that may be thrown by
	 *  individual implementations. For example, a SCOP identifier may require
	 *  that the domain definitions be available for download.
	 */
	public StructureIdentifier getBaseIdentifier() throws StructureException {
		if( base == null ) {

			switch(mySource) {
			case CATH:
				base = CathFactory.getCathDatabase().getDescriptionByCathId(getIdentifier());
				break;
			case ECOD:
				try {
					base = EcodFactory.getEcodDatabase().getDomainsById(name);
				} catch (IOException e) {
					throw new StructureException("Unable to get ECOD domain "+name,e);
				}
				break;
			case SCOP:
				// Fuzzy matching of the domain name to the current default factory
				base = guessScopDomain(getIdentifier(),ScopFactory.getSCOP());
				if(base == null) {
					// Guessing didn't work, so just use the PDBID and Chain from name
					// Guess that '_' means 'whole structure'
					if (chainName.equals("_")) {
						base = new SubstructureIdentifier(pdbId);
					} else {
						base = new SubstructureIdentifier(pdbId,ResidueRange.parseMultiple(chainName));
					}
					logger.error("Unable to find {}, so using {}",name,base);
				}
				break;
			case FILE:
				try {
					String[] prefix = name.split(":", 2);
					String filename;
					if(prefix.length > 1) {
						filename = prefix[1];
					} else {
						filename = name;
					}
					filename = FileDownloadUtils.expandUserHome(filename);
					base = new URLIdentifier(new File(filename).toURI().toURL());
				} catch (MalformedURLException e) {
					// Should never happen
					throw new StructureException("Unable to get URL for file: "+name,e);
				}
				break;
			case URL:
				try {
					base = new URLIdentifier(name);
				} catch (MalformedURLException e) {
					throw new StructureException("Invalid URL: "+name,e);
				}
				break;
			case PDP:
				try {
					PDPProvider provider = new RemotePDPProvider(false);
					base = provider.getPDPDomain(name);
				} catch (IOException e) {
					throw new StructureException("Unable to fetch PDP domain "+name, e);
				}
				break;
			case BIO:
				base = new BioAssemblyIdentifier(name);
				break;
			case PDB:
				base = new SubstructureIdentifier(getIdentifier());
				break;
			default:
				throw new IllegalStateException("Unimplemented source: "+mySource);
			}
		}
		return base;
	}

	@Override
	public SubstructureIdentifier toCanonical() throws StructureException {
		return getBaseIdentifier().toCanonical();
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return getBaseIdentifier().reduce(input);
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
	IOException {
		return getBaseIdentifier().loadStructure(cache);
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

	/**
	 * <p>
	 * Guess a scop domain. If an exact match is found, return that.
	 *
	 * <p>
	 * Otherwise, return the first scop domain found for the specified protein such that
	 * <ul>
	 * <li>The chains match, or one of the chains is '_' or '.'.
	 * <li>The domains match, or one of the domains is '_'.
	 * </ul>
	 *
	 * In some cases there may be several valid matches. In this case a warning
	 * will be logged.
	 *
	 * @param name SCOP domain name, or a guess thereof
	 * @param scopDB SCOP domain provider
	 * @return The best match for name among the domains of scopDB, or null if none match.
	 */
	public static ScopDomain guessScopDomain(String name, ScopDatabase scopDB) {
		List<ScopDomain> matches = new LinkedList<ScopDomain>();

		// Try exact match first
		ScopDomain domain = scopDB.getDomainByScopID(name);
		if (domain != null) {
			return domain;
		}

		// Didn't work. Guess it!
		logger.warn("Warning, could not find SCOP domain: " + name);

		Matcher scopMatch = scopPattern.matcher(name);
		if (scopMatch.matches()) {
			String pdbID = scopMatch.group(1);
			String chainName = scopMatch.group(2);
			String domainID = scopMatch.group(3);

			for (ScopDomain potentialSCOP : scopDB.getDomainsForPDB(pdbID)) {
				Matcher potMatch = scopPattern.matcher(potentialSCOP.getScopId());
				if (potMatch.matches()) {
					if (chainName.equals(potMatch.group(2)) || chainName.equals("_") || chainName.equals(".")
							|| potMatch.group(2).equals("_") || potMatch.group(2).equals(".")) {
						if (domainID.equals(potMatch.group(3)) || domainID.equals("_") || potMatch.group(3).equals("_")) {
							// Match, or near match
							matches.add(potentialSCOP);
						}
					}
				}
			}
		}

		Iterator<ScopDomain> match = matches.iterator();
		if (match.hasNext()) {
			ScopDomain bestMatch = match.next();
			if(logger.isWarnEnabled()) {
				StringBuilder warnMsg = new StringBuilder();
				warnMsg.append("Trying domain " + bestMatch.getScopId() + ".");
				if (match.hasNext()) {
					warnMsg.append(" Other possibilities: ");
					while (match.hasNext()) {
						warnMsg.append(match.next().getScopId()).append(" ");
					}
				}
				warnMsg.append(System.getProperty("line.separator"));
				logger.warn(warnMsg.toString());
			}
			return bestMatch;
		} else {
			return null;
		}
	}



}
