/*
 * BioJava development code
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
 * Author: Daniel Asarnow
 * Date:   2012-6-23
 */

package org.biojava.bio.structure.cath;

import java.util.List;

/** General API for interacting with CATH.
 *
 * @author Daniel Asarnow
 */
public interface CathDatabase {

    /** Return the CATH release version.
     *
     * @return CATH version
     */
    public String getCathVersion();

    /** Return the CathNode for a node ID.
     *
     * @param nodeId
     * @return CATH node
     */
    public CathNode getCathNode(String nodeId);

    /** Return list of CATH descriptions for node representatives at a CATH category (e.g. "T").
     *
     * @param category
     * @return CATH descriptions
     */
    public List<CathDomain> getByCategory(CathCategory category);

    /** Return list of CATH descriptions whose CATH codes (e.g. 1.4.6.10) start with the query.
     * This is currently redundant with getDescriptionsByNodeId.
     *
     * @param query
     * @return CATH descriptions
     */
    public List<CathDomain> filterByCathCode(String query);

    /** Return the CATH sub-tree for a particular domain.
     *
     * @param domain
     * @return CATH sub-tree
     */
    public List<CathNode> getTree(CathDomain domain);

    /** Return list of CATH domains whose node name (e.g. Orthogonal Bundle) starts with the query.
     *
     * @param query
     * @return CATH domains
     */
    public List<CathDomain> filterByNodeName(String query);

    /** Return list of CATH descriptions whose descriptions (name field) starts with the query.
     *
     * @param query
     * @return CATH descriptions
     */
    public List<CathDomain> filterByDescription(String query);

    /** Return CATH description for node representative by node ID.
     *
     * @param nodeId
     * @return CATH description
     */
    public CathDomain getDescriptionByNodeId(String nodeId);

    /** Return all CATH domains for a PDB ID.
     *
     * @param pdbId
     * @return CATH domains
     */
    public List<CathDomain> getDomainsForPdb(String pdbId);

    /** Return CATH domain for CATH domain ID.
     *
     * @param cathId
     * @return CATH domain
     */
    public CathDomain getDomainByCathId(String cathId);

    /** Return CATH description for CATH domain ID.
     *
     * @param cathId
     * @return
     */
    public CathDomain getDescriptionByCathId(String cathId);

    /** Return all CATH domains for a particular CATH node.
     *
     * @param nodeId
     * @return
     */
    public List<CathDomain> getDomainsByNodeId(String nodeId);

    public List<CathFragment> getFragmentsByPdbId(String pdbId);
}
