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
 */
package org.biojava.nbio.structure.ecod;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * An EcodDomain contains all the information of the ECOD database: id, 
 * classification groups (from higher to lower in the tree: X,H,T,F), PDB code,
 * chain, residue ranges and status (manual or automatic classification).
 * <p>
 * For detailed explanation about the ECOD information see the original article
 * at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4256011.
 * <pre>
 * Cheng H, Schaeffer RD, Liao Y, et al. 
 * ECOD: An Evolutionary Classification of Protein Domains. 
 * Elofsson A, ed. PLoS Computational Biology. 2014;10(12):e1003926.
 * </pre>
 * 
 * @author Spencer Bliven
 *
 */
public class EcodDomain implements Serializable, Cloneable, StructureIdentifier {

	/*
Column 1: ECOD uid - internal domain unique identifier
Column 2: ECOD domain id - domain identifier
Column 3: ECOD representative status - manual (curated) or automated nonrep
Column 4: ECOD hierachy identifier - [X-group].[H-group].[T-group].[F-group]
Column 5: PDB identifier
Column 6: Chain identifier (note: case-sensitive)
Column 7: PDB residue number range
Column 8: seq_id number range (based on internal PDB indices)
Column 9: Architecture name
Column 10: X-group name
Column 11: H-group name
Column 12: T-group name
Column 13: F-group name (F_UNCLASSIFIED denotes that domain has not been assigned to an F-group)
Column 14: Domain assembly status (if domain is member of assembly, partners' ecod domain ids listed)
Column 15: Comma-separated value list of non-polymer entities within 4 A of at least one residue of domain


001502751       e4s1gA1 1.1.1   4s1g    A       A:68-251        beta barrels    cradle loop barrel      RIFT-related    acid protease   F_UNCLASSIFIED  NOT_DOMAIN_ASSEMBLY     NO_LIGANDS_4A

	 */


	private static final long serialVersionUID = -7760082165560332048L;


	private Long uid;
	private String domainId;
	private Boolean manual;
	private Integer xGroup;
	private Integer hGroup;
	private Integer tGroup;
	private Integer fGroup;
	private String pdbId;
	private String chainId;
	private String range;
	private String seqIdRange;
	private String architectureName;
	private String xGroupName;
	private String hGroupName;
	private String tGroupName;
	private String fGroupName;
	private Long assemblyId; //for non-assemblies, matches the uid.
	private Set<String> ligands;

	/** Default constructor with all null properties */
	public EcodDomain() {}

	public EcodDomain(Long uid, String domainId, Boolean manual,
			Integer xGroup, Integer hGroup, Integer tGroup, Integer fGroup, String pdbId,
			String chainId, String range, String architectureName,
			String xGroupName, String hGroupName, String tGroupName,
			String fGroupName, Long assemblyId, Set<String> ligands) {
		this(uid, domainId, manual,
				xGroup, hGroup, tGroup, fGroup, pdbId,
				chainId, range, null, architectureName,
				xGroupName, hGroupName, tGroupName,
				fGroupName, assemblyId, ligands);
	}
	public EcodDomain(Long uid, String domainId, Boolean manual,
				Integer xGroup, Integer hGroup, Integer tGroup, Integer fGroup, String pdbId,
				String chainId, String range, String seqId, String architectureName,
				String xGroupName, String hGroupName, String tGroupName,
				String fGroupName, Long assemblyId, Set<String> ligands) {
		this.uid = uid;
		this.domainId = domainId;
		this.manual = manual;
		this.xGroup = xGroup;
		this.hGroup = hGroup;
		this.tGroup = tGroup;
		this.fGroup = fGroup;
		this.pdbId = pdbId;
		this.chainId = chainId;
		this.range = range;
		this.seqIdRange = seqId;
		this.architectureName = architectureName;
		this.xGroupName = xGroupName;
		this.hGroupName = hGroupName;
		this.tGroupName = tGroupName;
		this.fGroupName = fGroupName;
		this.assemblyId = assemblyId;
		this.ligands = ligands;
	}
	public EcodDomain(String domainId) {
		this.domainId = domainId;
	}
	public EcodDomain(EcodDomain o) {
		this.uid = o.uid;
		this.domainId = o.domainId;
		this.manual = o.manual;
		this.xGroup = o.xGroup;
		this.hGroup = o.hGroup;
		this.tGroup = o.tGroup;
		this.fGroup = o.fGroup;
		this.pdbId = o.pdbId;
		this.chainId = o.chainId;
		this.range = o.range;
		this.seqIdRange = o.seqIdRange;
		this.architectureName = o.architectureName;
		this.xGroupName = o.xGroupName;
		this.hGroupName = o.hGroupName;
		this.tGroupName = o.tGroupName;
		this.fGroupName = o.fGroupName;
		this.assemblyId = o.assemblyId;
		this.ligands = new HashSet<String>(o.ligands);
	}



	@Override
	protected Object clone() throws CloneNotSupportedException {
		return new EcodDomain(this);
	}

	public Long getUid() {
		return uid;
	}
	public void setUid(Long uid) {
		this.uid = uid;
	}
	public String getDomainId() {
		return domainId;
	}
	public void setDomainId(String domainId) {
		this.domainId = domainId;
	}
	public Boolean getManual() {
		return manual;
	}
	public void setManual(Boolean manual) {
		this.manual = manual;
	}
	public Integer getXGroup() {
		return xGroup;
	}
	public void setXGroup(Integer xGroup) {
		this.xGroup = xGroup;
	}
	public Integer getHGroup() {
		return hGroup;
	}
	public void setHGroup(Integer hGroup) {
		this.hGroup = hGroup;
	}
	public Integer getTGroup() {
		return tGroup;
	}
	public void setTGroup(Integer tGroup) {
		this.tGroup = tGroup;
	}
	public Integer getFGroup() {
		return fGroup;
	}
	public void setFGroup(Integer fGroup) {
		this.fGroup = fGroup;
	}
	public String getPdbId() {
		return pdbId;
	}
	public void setPdbId(String pdbId) {
		this.pdbId = pdbId;
	}
	public String getChainId() {
		return chainId;
	}
	public void setChainId(String chainId) {
		this.chainId = chainId;
	}
	/**
	 * Get the range of this domain, in PDB residue numbers (mmCif's
	 * _pdbx_poly_seq_scheme.pdb_seq_num and pdb_ins_code).
	 * @return The chain and residue range, e.g. "A:1-100"
	 */
	public String getRange() {
		return range;
	}
	public void setRange(String range) {
		this.range = range;
	}
	/**
	 * Get the range of this domain, in 1-based residue indices (mmCif's
	 * _pdbx_poly_seq_scheme.seq_id)
	 *
	 * Note that {@link #getRange()} is used when constructing the domain.
	 * @return The chain and residue range, e.g. "A:1-100"
	 */
	public String getSeqIdRange() {
		return seqIdRange;
	}
	public void setSeqIdRange(String seqIdRange) {
		this.seqIdRange = seqIdRange;
	}
	public String getArchitectureName() {
		return architectureName;
	}
	public void setArchitectureName(String architectureName) {
		this.architectureName = architectureName;
	}
	public String getXGroupName() {
		return xGroupName;
	}
	public void setXGroupName(String xGroupName) {
		this.xGroupName = xGroupName;
	}
	public String getHGroupName() {
		return hGroupName;
	}
	public void setHGroupName(String hGroupName) {
		this.hGroupName = hGroupName;
	}
	public String getTGroupName() {
		return tGroupName;
	}
	public void setGroupName(String tGroupName) {
		this.tGroupName = tGroupName;
	}
	public String getFGroupName() {
		return fGroupName;
	}
	public void setFGroupName(String fGroupName) {
		this.fGroupName = fGroupName;
	}
	/**
	 * @return The assembly ID, or the DomainId if not in an assembly, or null if unknown.
	 */
	public Long getAssemblyId() {
		return assemblyId;
	}
	public void setAssemblyId(Long assemblyId) {
		this.assemblyId = assemblyId;
	}
	public Set<String> getLigands() {
		return ligands;
	}
	public void setLigands(Set<String> ligands) {
		this.ligands = ligands;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "EcodDomain [uid=" + uid + ", domainId=" + domainId
				+ ", manual=" + manual + ", xGroup=" + xGroup + ", hGroup="
				+ hGroup + ", tGroup=" + tGroup + ", fGroup="+ fGroup + ", pdbId=" + pdbId
				+ ", chainName=" + chainId + ", range=" + range
				+ ", architectureName=" + architectureName + ", xGroupName="
				+ xGroupName + ", hGroupName=" + hGroupName + ", tGroupName="
				+ tGroupName + ", fGroupName=" + fGroupName + ", assemblyId="
				+ assemblyId + ", ligands=" + ligands + "]";
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime
				* result
				+ ((architectureName == null) ? 0 : architectureName.hashCode());
		result = prime * result + ((chainId == null) ? 0 : chainId.hashCode());
		result = prime * result
				+ ((domainId == null) ? 0 : domainId.hashCode());
		result = prime * result
				+ ((fGroupName == null) ? 0 : fGroupName.hashCode());
		result = prime * result + ((fGroup == null) ? 0 : fGroup.hashCode());
		result = prime * result + ((hGroup == null) ? 0 : hGroup.hashCode());
		result = prime * result
				+ ((hGroupName == null) ? 0 : hGroupName.hashCode());
		result = prime * result
				+ ((assemblyId == null) ? 0 : assemblyId.hashCode());
		result = prime * result + ((ligands == null) ? 0 : ligands.hashCode());
		result = prime * result + ((manual == null) ? 0 : manual.hashCode());
		result = prime * result + ((pdbId == null) ? 0 : pdbId.hashCode());
		result = prime * result + ((range == null) ? 0 : range.hashCode());
		result = prime * result + ((tGroup == null) ? 0 : tGroup.hashCode());
		result = prime * result
				+ ((tGroupName == null) ? 0 : tGroupName.hashCode());
		result = prime * result + ((uid == null) ? 0 : uid.hashCode());
		result = prime * result + ((xGroup == null) ? 0 : xGroup.hashCode());
		result = prime * result
				+ ((xGroupName == null) ? 0 : xGroupName.hashCode());
		return result;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		EcodDomain other = (EcodDomain) obj;
		if (architectureName == null) {
			if (other.architectureName != null)
				return false;
		} else if (!architectureName.equals(other.architectureName))
			return false;
		if (chainId == null) {
			if (other.chainId != null)
				return false;
		} else if (!chainId.equals(other.chainId))
			return false;
		if (domainId == null) {
			if (other.domainId != null)
				return false;
		} else if (!domainId.equals(other.domainId))
			return false;
		if (fGroupName == null) {
			if (other.fGroupName != null)
				return false;
		} else if (!fGroupName.equals(other.fGroupName))
			return false;
		if (fGroup == null) {
			if (other.fGroup != null)
				return false;
		} else if (!fGroup.equals(other.fGroup))
			return false;
		if (hGroup == null) {
			if (other.hGroup != null)
				return false;
		} else if (!hGroup.equals(other.hGroup))
			return false;
		if (hGroupName == null) {
			if (other.hGroupName != null)
				return false;
		} else if (!hGroupName.equals(other.hGroupName))
			return false;
		if (assemblyId == null) {
			if (other.assemblyId != null)
				return false;
		} else if (!assemblyId.equals(other.assemblyId))
			return false;
		if (ligands == null) {
			if (other.ligands != null)
				return false;
		} else if (!ligands.equals(other.ligands))
			return false;
		if (manual == null) {
			if (other.manual != null)
				return false;
		} else if (!manual.equals(other.manual))
			return false;
		if (pdbId == null) {
			if (other.pdbId != null)
				return false;
		} else if (!pdbId.equals(other.pdbId))
			return false;
		if (range == null) {
			if (other.range != null)
				return false;
		} else if (!range.equals(other.range))
			return false;
		if (tGroup == null) {
			if (other.tGroup != null)
				return false;
		} else if (!tGroup.equals(other.tGroup))
			return false;
		if (tGroupName == null) {
			if (other.tGroupName != null)
				return false;
		} else if (!tGroupName.equals(other.tGroupName))
			return false;
		if (uid == null) {
			if (other.uid != null)
				return false;
		} else if (!uid.equals(other.uid))
			return false;
		if (xGroup == null) {
			if (other.xGroup != null)
				return false;
		} else if (!xGroup.equals(other.xGroup))
			return false;
		if (xGroupName == null) {
			if (other.xGroupName != null)
				return false;
		} else if (!xGroupName.equals(other.xGroupName))
			return false;
		return true;
	}

	@Override
	public String getIdentifier() {
		return getDomainId();
	}

	public List<ResidueRange> getResidueRanges() {
		return ResidueRange.parseMultiple(range);
	}

	@Override
	public SubstructureIdentifier toCanonical() {
		return new SubstructureIdentifier(getPdbId(), ResidueRange.parseMultiple(getRange()));
	}

	@Override
	public Structure reduce(Structure input) throws StructureException {
		return toCanonical().reduce(input);
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return cache.getStructureForPdbId(pdbId);
	}

}
