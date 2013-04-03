package org.biojava.bio.structure.scop;

import java.util.List;

/** General API how to interact with SCOP
 * 
 * @author Andreas Prlic
 * @since 3.0.2
 */
public interface ScopDatabase {

	/** Get all records of a particular classification.
	 * 
	 * @param category e.g. "superfamily"
	 * @return all records of this type
	 */
	public abstract List<ScopDescription> getByCategory(ScopCategory category);

	/** Get all scop descriptions that start with a classifcation ID, e.g. b.1.18
	 * 
	 * @param query
	 * @return list of scop descriptions
	 */
	public abstract List<ScopDescription> filterByClassificationId(String query);

	/** get the SCOP sub-tree for a particular domain. 
	 * 
	 * @param domain
	 * @return list of ScopNodes providing the path to this domain
	 */
	public abstract List<ScopNode> getTree(ScopDomain domain);

	
	/** search through SCOP and filter based on domain name
	 * 
	 * @param query a (part) of a name
	 * @return list of matchin ScopDomains
	 */
	public abstract List<ScopDomain> filterByDomainName(String query);

	/** Get all scop descriptions that start with a certain name. e.g. Globin
	 * 
	 * @param query
	 * @return list of scop descriptions
	 */
	public abstract List<ScopDescription> filterByDescription(String query);

	/** Return the SCOP description for a node in the hierarchy by its "sunid" id.
	 * 
	 * @param sunid
	 * @return a ScopDescription object
	 */
	public abstract ScopDescription getScopDescriptionBySunid(int sunid);

	/** Get a list of ScopDomains that have been assigned to a PDB ID
	 * 
	 * @param pdbId the PDB entry
	 * @return a list of ScopDomains
	 */
	public abstract List<ScopDomain> getDomainsForPDB(String pdbId);

	/** get a ScopDomain by its SCOP ID (warning, they are not stable between releases!)
	 * 
	 *
	 * @param scopId e.g. d2bq6a1 
	 * @return a ScopDomain or null if no domain with the particular ID could be found
	 */
	public abstract ScopDomain getDomainByScopID(String scopId);

	/** Access a particular ScopNode. The scopNode then allows to traverse through the scop hierarchy...
	 * 
	 * @param sunid the scop unique id
	 * @return a ScopNode that matches this sunid
	 */
	public abstract ScopNode getScopNode(int sunid);
	
	/** Returns the SCOP version
	 * 
	 * @return version of SCOP
	 */

	public abstract String getScopVersion();

	/** Get a SCOP domain by its sunid
	 * 
	 * @param sunid the scop unique id
	 * @return a list of scop domains that match this sunid
	 */
	public abstract List<ScopDomain> getScopDomainsBySunid(Integer sunid);

}