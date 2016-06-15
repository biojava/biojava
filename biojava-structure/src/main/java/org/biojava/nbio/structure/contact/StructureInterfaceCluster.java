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
package org.biojava.nbio.structure.contact;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class StructureInterfaceCluster implements Serializable {

	private static final long serialVersionUID = 1L;


	private int id;

	private List<StructureInterface> members;

	/**
	 * The average similarity score between all pairs of members in the cluster.
	 */
	private double averageScore;


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

	/**
	 * Returns the average similarity score between all pairs of members in the cluster
	 * @return
	 */
	public double getAverageScore() {
		return averageScore;
	}

	public void setAverageScore(double averageScore) {
		this.averageScore = averageScore;
	}
}
