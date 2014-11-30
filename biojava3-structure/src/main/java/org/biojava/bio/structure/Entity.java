package org.biojava.bio.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * An entity, i.e. a group of sequence identical chains 
 * related by non-crystallographic symmetry (NCS).
 * 
 * @author duarte_j
 *
 */
public class Entity implements Serializable {

	private static final long serialVersionUID = 1L;

	/**
	 * The representative chain
	 */
	private Chain representative;
	
	/**
	 * The member chains including the representative
	 */
	private List<Chain> members;
	
	public Entity() {
		members = new ArrayList<Chain>();
	}
	
	/**
	 * Construct a new Entity given representative and members
	 * @param representative
	 * @param members
	 * @throws IllegalArgumentException if the members don't contain the representative
	 */
	public Entity(Chain representative, List<Chain> members) {
		if (!members.contains(representative))
			throw new IllegalArgumentException("The member chains given don't contain the representative chain");
		
		this.representative = representative;		
		this.members = members;
	}

	/**
	 * Get the representative chain
	 * @return
	 */
	public Chain getRepresentative() {
		return representative;
	}
	
	public void setRepresentative(Chain representative) {
		this.representative = representative;
	}
	
	/**
	 * Get the member chains including the representative
	 * @return
	 */
	public List<Chain> getMembers() {
		return members;
	}
	
	public void setMembers(List<Chain> members) {
		this.members = members;
	}
	
	public boolean addMember(Chain c) {
		return members.add(c);		
	}
}
