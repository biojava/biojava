/**
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
 * Created on Nov 27, 2012
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.symmetry.analysis;


import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.symmetry.core.AxisAligner;
import org.biojava.bio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.bio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.bio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.bio.structure.symmetry.core.RotationGroup;
import org.biojava.bio.structure.symmetry.core.Subunits;
import org.biojava.bio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;


public class CalcBioAssemblySymmetry {
	private Structure bioAssembly;
	private QuatSymmetryParameters parameters;
	private QuatSymmetryResults results;

	private JmolSymmetryScriptGenerator scriptGenerator;
	
	static public String version;
	static public String build; 
	static {
		try {
			ResourceManager about = ResourceManager.getResourceManager("about");

			version = about.getString("project_version");
			build   = about.getString("build");

		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public CalcBioAssemblySymmetry(Structure bioAssembly, QuatSymmetryParameters parameters){
		this.bioAssembly = bioAssembly;
		this.parameters = parameters;
	}
	
	public QuatSymmetryParameters getParameters(){
		return parameters;
	}

	public QuatSymmetryDetector orient(){
		QuatSymmetryDetector detector = new QuatSymmetryDetector(bioAssembly, parameters);
		String defaultColoring = "";

		if (detector.hasProteinSubunits()) {	
			for (QuatSymmetryResults globalSymmetry: detector.getGlobalSymmetry()) {

				String postFix = "g";
				AxisAligner aligner = AxisAligner.getInstance(globalSymmetry);
				JmolSymmetryScriptGenerator generator = JmolSymmetryScriptGenerator.getInstance(aligner, postFix);
				generator.setOnTheFly(parameters.isOnTheFly());
				// save the preferred result
				if (globalSymmetry.isPreferredResult()) {
					this.results = globalSymmetry;
					this.scriptGenerator = generator;
				}
					
				// get default color by symmetry for all subunits (for asymmetric cases only)
				if (globalSymmetry.getSymmetry().equals("C1") && ! globalSymmetry.getSubunits().isPseudoStoichiometric()) {
					defaultColoring = generator.colorBySymmetry();
				}
				
				if (parameters.isVerbose()) {
					System.out.println("Global symmetry: ");
					System.out.println(globalSymmetry);
					
//					System.out.println();
//					System.out.println(generator.getDefaultOrientation());
//					System.out.println(generator.getZoom());
//					System.out.println(generator.drawPolyhedron());
					System.out.println(generator.drawAxes());
					System.out.println(generator.colorBySubunit());
					System.out.println(generator.colorBySequenceCluster());
					System.out.println(generator.colorBySymmetry());
				}
				
			}

			for (List<QuatSymmetryResults> localSymmetries: detector.getLocalSymmetries()) {	
				int count = 0;

				for (QuatSymmetryResults localSymmetry: localSymmetries) {
					// create a unique postFix for each local symmetry to be used in the Jmol script
					// to avoid naming conflicts when more than one symmetry is displayed at one time.
					String postFix = "l" + count;
					AxisAligner aligner = AxisAligner.getInstance(localSymmetry);
					JmolSymmetryScriptGenerator generator = JmolSymmetryScriptGenerator.getInstance(aligner, postFix);		
					generator.setOnTheFly(parameters.isOnTheFly());
					// sets color by symmetry for all subunits. This is
					// required for local symmetry, to ensure all subunits are colored.
					generator.setDefaultColoring(defaultColoring);
					// save preferred result
					if (localSymmetry.isPreferredResult()) {
						this.results = localSymmetry;
						this.scriptGenerator = generator;
					}

					if (parameters.isVerbose()) {
						System.out.println("Local symmetry: ");
						System.out.println(localSymmetry);
						
//						System.out.println();
//						System.out.println(generator.getDefaultOrientation());
//						System.out.println(generator.getZoom());
//						System.out.println(generator.drawPolyhedron());
//						System.out.println(generator.drawAxes());
//						System.out.println(generator.colorBySubunit());
//						System.out.println(generator.colorBySequenceCluster());
//						System.out.println(generator.colorBySymmetry());
					}
				}
			}
		} else {
			System.out.println("No protein chains found for " + bioAssembly.getPDBCode() );
		}
		return detector;
	}

	/** Only works if this is not helical symmetry. Deprecated, use getSymmetry instead */
	@Deprecated
	public RotationGroup getRotationGroup() {
		return results.getRotationGroup();
	}

	
	public Subunits getSubunits() {
		return results.getSubunits();
	}
	
	/** String representation of Symmetry, e.g. C1, C2, H 
	 * 
	 * @return
	 */
	public String getSymmetry() {
		return results.getSymmetry();
	}

	public JmolSymmetryScriptGenerator getScriptGenerator() {
		return scriptGenerator;
	}
}
