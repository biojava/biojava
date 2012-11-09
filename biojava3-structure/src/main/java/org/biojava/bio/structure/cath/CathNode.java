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

import java.io.Serializable;

/** Represents a node in the CATH hierarchy.
 *
 * @author Daniel Asarnow
 */
public class CathNode implements Serializable{

    public static final long serialVersionUID = 1L;

    /**
     * The CATH code of the node, e.g. "1.10.8.10".
     */
    String nodeId;

    /**
     * The CATH code of the parent, e.g. "1.10.8". Calculated during parsing.
     */
    String parentId;

    /**
     * The representative domain for this node.
     */
    String representative;

    /**
     * A name or description.
     */
    String description;

    /**
     * This node's level within the hierarchy. One of CATH, not CATHSOLID.
     */
    CathCategory category;

    public String getNodeId() {
        return nodeId;
    }

    public void setNodeId(String nodeId) {
        this.nodeId = nodeId;
        this.category = CathCategory.fromCathCode(nodeId);
    }

    public String getParentId() {
        return parentId;
    }

    public void setParentId(String parentId) {
        this.parentId = parentId;
    }

    public String getRepresentative() {
        return representative;
    }

    public void setRepresentative(String representative) {
        this.representative = representative;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public CathCategory getCategory() {
        return category;
    }

    @Override
    public String toString() {
        return ""; //TODO implement
    }
}
