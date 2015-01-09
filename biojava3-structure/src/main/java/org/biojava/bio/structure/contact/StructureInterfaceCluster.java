package org.biojava.bio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class StructureInterfaceCluster implements Serializable {

	private static final long serialVersionUID = 1L;
	
	
	private int id;
	
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
	
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	/**
	 * Return the average buried surface area for this interface cluster
	 * @return
	 */
	public double getTotalArea() {
		double area = 0;
		for (StructureInterface interf:members) {
			area+=interf.getTotalArea();
		}
		return area/members.size();
	}
}
