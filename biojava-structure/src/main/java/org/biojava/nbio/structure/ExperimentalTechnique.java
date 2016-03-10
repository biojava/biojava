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
package org.biojava.nbio.structure;

import java.util.HashMap;
import java.util.Set;

/**
 *
 * An enum to represent the experimental technique of a PDB structure
 *
 * @author duarte_j
 *
 */
public enum ExperimentalTechnique {


	 XRAY_DIFFRACTION			("X-RAY DIFFRACTION", 		true,	false),

	 SOLUTION_NMR				("SOLUTION NMR",			false,	true),
	 SOLID_STATE_NMR			("SOLID-STATE NMR",			false,	true),

	 ELECTRON_MICROSCOPY		("ELECTRON MICROSCOPY",		false,	false),
	 ELECTRON_CRYSTALLOGRAPHY	("ELECTRON CRYSTALLOGRAPHY",true,	false),

	 FIBER_DIFFRACTION			("FIBER DIFFRACTION",		false,	false),

	 NEUTRON_DIFFRACTION		("NEUTRON DIFFRACTION",		true,	false),

	 SOLUTION_SCATTERING		("SOLUTION SCATTERING",		false,	false),


	 // from here, not in "official list" in pdb.org advanced search: they call these "OTHER"
	 POWDER_DIFFRACTION			("POWDER DIFFRACTION",		true,	false),

	 FLUORESCENCE_TRANSFER		("FLUORESCENCE TRANSFER",	false,	false),

	 INFRARED_SPECTROSCOPY		("INFRARED SPECTROSCOPY",	false,	false);


	 private static final HashMap<String, ExperimentalTechnique> expTechStr2Value = initExpTechStr2Value();


	 private String name;
	 private boolean isCrystallographic;
	 private boolean isNmr;

	 private ExperimentalTechnique(String name, boolean isXtallographic, boolean isNmr) {
		 this.name = name;
		 this.isCrystallographic = isXtallographic;
		 this.isNmr = isNmr;
	 }


	 private static HashMap<String, ExperimentalTechnique> initExpTechStr2Value() {
		HashMap<String, ExperimentalTechnique> expTechStr2Value = new HashMap<String, ExperimentalTechnique>();
		for(ExperimentalTechnique exp:ExperimentalTechnique.values()) {
			expTechStr2Value.put(exp.getName(), exp);
		}
		return expTechStr2Value;
	 }

	 public String getName() {
		 return name;
	 }

	 public boolean isCrystallographic() {
		 return isCrystallographic;
	 }

	 public boolean isNmr() {
		 return isNmr;
	 }

	 /**
	  * Returns the ExpTechnique given an experimental technique name as used in the PDB,
	  * e.g. "X-RAY DIFFRACTION" returns {@link ExperimentalTechnique#XRAY_DIFFRACTION}
	  * @param expTechniqueName the ExpTechnique value or null if string doesn't match one of the known PDB experimental strings
	  * @return
	  */
	 public static ExperimentalTechnique getByName(String expTechniqueName) {
		 return expTechStr2Value.get(expTechniqueName);
	 }

	 /**
	  * Given a Set of ExperimentalTechniques returns true if at least one is crystallographic
	  * @return true if at least 1 of the techniques is crystallographic, false if
	  * none of the techniques are crystallographic
	  * @throws NullPointerException if input is null
	  */
	 public static boolean isCrystallographic(Set<ExperimentalTechnique> techniques) {

		 for (ExperimentalTechnique et:techniques) {
			 if (et.isCrystallographic()) return true;
		 }

		 return false;

	 }

	 /**
	  * Given a Set of ExperimentalTechniques returns true if at least one is NMR
	  * @return true if at least 1 of the techniques is NMR, false if
	  * none of the techniques are NMR
	  * @throws NullPointerException if input is null
	  */
	 public static boolean isNmr(Set<ExperimentalTechnique> techniques) {

		 for (ExperimentalTechnique et:techniques) {
			 if (et.isNmr()) return true;
		 }

		 return false;
	 }

	 @Override
	public String toString() {
		 return getName();
	 }
}
