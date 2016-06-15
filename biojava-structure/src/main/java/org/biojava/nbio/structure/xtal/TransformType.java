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

public enum TransformType {

	//              id fold  screw infinite  shortName
	AU				(0,   1, false, false,   "AU"),
	XTALTRANSL		(1,   1, false, true,    "XT"),  // translation
	CELLTRANSL		(2,   1, false, true,    "FT"),  // fractional translation

	TWOFOLD			(3,   2, false, false,   "2" ),
	TWOFOLDSCREW	(4,   2, true , true,    "2S"),

	THREEFOLD		(5,   3, false, false,   "3" ),
	THREEFOLDSCREW	(6,   3, true,  true,    "3S"),

	FOURFOLD		(7,   4, false, false,   "4" ),
	FOURFOLDSCREW	(8,   4, true,  true,    "4S"),

	SIXFOLD       	(9,   6, false, false,   "6" ),
	SIXFOLDSCREW  	(10,  6, true,  true,    "6S"),

	ONEBAR          (11, -1, false, false,   "-1"),

	TWOBAR          (12, -2, false, false,   "-2"),
	GLIDE           (13, -2, true,  false,   "GL"),

	THREEBAR        (14, -3, false, false,   "-3"),

	FOURBAR         (15, -4, false, false,   "-4"),

	SIXBAR          (16, -6, false, false,   "-6");



	private int id;
	private int foldType;
	private boolean isScrew;
	private boolean isInfinite;
	private String shortName;

	private TransformType(int id, int foldType, boolean isScrew, boolean isInfinite, String shortName) {
		this.id = id;
		this.foldType = foldType;
		this.isScrew = isScrew;
		this.isInfinite = isInfinite;
		this.shortName = shortName;
	}

	public int getId() {
		return id;
	}

	public int getFoldType() {
		return foldType;
	}

	/**
	 * Tells whether the transform type is a screw or glide plane
	 * @return
	 */
	public boolean isScrew() {
		return isScrew;
	}

	/**
	 * Tells whether the transform type produces infinite assemblies
	 * if interface happens between identical chains
	 * @return
	 */
	public boolean isInfinite() {
		return isInfinite;
	}

	public String getShortName() {
		return shortName;
	}
}
