package org.biojava.bio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class StructureInterfaceCluster implements Serializable {

	private static final long serialVersionUID = 1L;
	
	
	private List<StructureInterface> members;
	
	
	
	public StructureInterfaceCluster() {
		this.members = new ArrayList<StructureInterface>();
	}
	
	public List<StructureInterface> getMembers() {
		return members;
	}
	
	public void setMembers(List<StructureInterface> members) {
		this.members = members;
	}

	public boolean addMember(StructureInterface interf) {
		return this.members.add(interf);
	}
}
