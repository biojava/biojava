package org.biojava.nbio.structure.domain;

import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class EcodDomain implements Serializable, Cloneable {

	/*
Column 1: ECOD uid - internal domain unique identifier
Column 2: ECOD domain id - domain identifier
Column 3: ECOD representative status - manual (curated) or automated nonrep
Column 4: ECOD hierachy identifier - [X-group].[H-group].{T-group]
Column 5: PDB identifier
Column 6: Chain identifier (note: case-sensitive)
Column 7: PDB residue number range
Column 8: Architecture name
Column 9: X-group name
Column 10: H-group name
Column 11: T-group name
Column 12: F-group name (F_UNCLASSIFIED denotes that domain has not been assigned to an F-group)
Column 13: Domain assembly status (if domain is member of assembly, partners' ecod domain ids listed)
Column 14: Comma-separated value list of non-polymer entities within 4 A of at least one residue of domain


001502751       e4s1gA1 1.1.1   4s1g    A       A:68-251        beta barrels    cradle loop barrel      RIFT-related    acid protease   F_UNCLASSIFIED  NOT_DOMAIN_ASSEMBLY     NO_LIGANDS_4A

	 */


	private static final long serialVersionUID = -7760082165560332048L;

	/** String for unclassified F-groups */
	public static final String F_UNCLASSFIED = "F_UNCLASSIFIED";

	private Long uid;
	private String domainId;
	private Boolean manual;
	private Integer xGroup;
	private Integer hGroup;
	private Integer tGroup;
	private String pdbId;
	private String chainId;
	private String range;
	private String architectureName;
	private String xGroupName;
	private String hGroupName;
	private String tGroupName;
	private String fGroupName;
	private Boolean isAssembly; // Maybe should be a list, according to description?
	private Set<String> ligands;
	
	/** Default constructor with all null properties */
	public EcodDomain() {}
	
	public EcodDomain(Long uid, String domainId, Boolean manual,
			Integer xGroup, Integer hGroup, Integer tGroup, String pdbId,
			String chainId, String range, String architectureName,
			String xGroupName, String hGroupName, String tGroupName,
			String fGroupName, Boolean isAssembly, Set<String> ligands) {
		this.uid = uid;
		this.domainId = domainId;
		this.manual = manual;
		this.xGroup = xGroup;
		this.hGroup = hGroup;
		this.tGroup = tGroup;
		this.pdbId = pdbId;
		this.chainId = chainId;
		this.range = range;
		this.architectureName = architectureName;
		this.xGroupName = xGroupName;
		this.hGroupName = hGroupName;
		this.tGroupName = tGroupName;
		this.fGroupName = fGroupName;
		this.isAssembly = isAssembly;
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
		this.pdbId = o.pdbId;
		this.chainId = o.chainId;
		this.range = o.range;
		this.architectureName = o.architectureName;
		this.xGroupName = o.xGroupName;
		this.hGroupName = o.hGroupName;
		this.tGroupName = o.tGroupName;
		this.fGroupName = o.fGroupName;
		this.isAssembly = o.isAssembly;
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
	public Integer getxGroup() {
		return xGroup;
	}
	public void setxGroup(Integer xGroup) {
		this.xGroup = xGroup;
	}
	public Integer gethGroup() {
		return hGroup;
	}
	public void sethGroup(Integer hGroup) {
		this.hGroup = hGroup;
	}
	public Integer gettGroup() {
		return tGroup;
	}
	public void settGroup(Integer tGroup) {
		this.tGroup = tGroup;
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
	public String getRange() {
		return range;
	}
	public void setRange(String range) {
		this.range = range;
	}
	public String getArchitectureName() {
		return architectureName;
	}
	public void setArchitectureName(String architectureName) {
		this.architectureName = architectureName;
	}
	public String getxGroupName() {
		return xGroupName;
	}
	public void setxGroupName(String xGroupName) {
		this.xGroupName = xGroupName;
	}
	public String gethGroupName() {
		return hGroupName;
	}
	public void sethGroupName(String hGroupName) {
		this.hGroupName = hGroupName;
	}
	public String gettGroupName() {
		return tGroupName;
	}
	public void settGroupName(String tGroupName) {
		this.tGroupName = tGroupName;
	}
	public String getfGroupName() {
		return fGroupName;
	}
	public void setfGroupName(String fGroupName) {
		this.fGroupName = fGroupName;
	}
	public Boolean getIsAssembly() {
		return isAssembly;
	}
	public void setIsAssembly(Boolean isAssembly) {
		this.isAssembly = isAssembly;
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
				+ hGroup + ", tGroup=" + tGroup + ", pdbId=" + pdbId
				+ ", chainId=" + chainId + ", range=" + range
				+ ", architectureName=" + architectureName + ", xGroupName="
				+ xGroupName + ", hGroupName=" + hGroupName + ", tGroupName="
				+ tGroupName + ", fGroupName=" + fGroupName + ", isAssembly="
				+ isAssembly + ", ligands=" + ligands + "]";
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
		result = prime * result + ((hGroup == null) ? 0 : hGroup.hashCode());
		result = prime * result
				+ ((hGroupName == null) ? 0 : hGroupName.hashCode());
		result = prime * result
				+ ((isAssembly == null) ? 0 : isAssembly.hashCode());
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
		if (isAssembly == null) {
			if (other.isAssembly != null)
				return false;
		} else if (!isAssembly.equals(other.isAssembly))
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

}
