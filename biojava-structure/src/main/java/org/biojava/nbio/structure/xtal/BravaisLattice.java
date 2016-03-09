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
package org.biojava.nbio.structure.xtal;

import java.util.HashMap;

/**
 * An enum to represent the 7 Bravais lattices
 *
 * @author duarte_j
 *
 */
public enum BravaisLattice {

	TRICLINIC    (1, "TRICLINIC",    new CrystalCell(1.00,1.25,1.50, 60,70,80)), // alpha,beta,gamma!=90
	MONOCLINIC   (2, "MONOCLINIC",   new CrystalCell(1.00,1.25,1.50, 90,60,90)), // beta!=90, alpha=gamma=90
	ORTHORHOMBIC (3, "ORTHORHOMBIC", new CrystalCell(1.00,1.25,1.50, 90,90,90)), // alpha=beta=gamma=90
	TETRAGONAL   (4, "TETRAGONAL",   new CrystalCell(1.00,1.00,1.25, 90,90,90)), // alpha=beta=gamma=90, a=b
	TRIGONAL     (5, "TRIGONAL",     new CrystalCell(1.00,1.00,1.25, 90,90,120)),// a=b!=c, alpha=beta=90, gamma=120
	HEXAGONAL    (6, "HEXAGONAL",    new CrystalCell(1.00,1.00,1.25, 90,90,120)),// a=b!=c, alpha=beta=90, gamma=120
	CUBIC        (7, "CUBIC",        new CrystalCell(1.00,1.00,1.00, 90,90,90)); // a=b=c, alpha=beta=gamma=90

	private static HashMap<String, BravaisLattice> name2bl = initname2bl();
	private String name;
	private int id;
	private CrystalCell exampleUnitCell;

	private BravaisLattice(int id, String name, CrystalCell exampleUnitCell) {
		this.name = name;
		this.id = id;
		this.exampleUnitCell = exampleUnitCell;
	}

	public String getName() {
		return name;
	}

	public int getId() {
		return id;
	}

	public CrystalCell getExampleUnitCell() {
		return exampleUnitCell;
	}

	private static HashMap<String,BravaisLattice> initname2bl(){
		HashMap<String,BravaisLattice> name2bl = new HashMap<String, BravaisLattice>();
		for (BravaisLattice bl:BravaisLattice.values()) {
			name2bl.put(bl.getName(), bl);
		}
		return name2bl;
	}

	public static BravaisLattice getByName(String blName) {
		return name2bl.get(blName);
	}
}
