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
package demo;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.symmetry.analysis.CalcBioAssemblySymmetry;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;

import java.io.IOException;
import java.util.List;

/**
 * Created by ap3 on 02/04/2015.
 */
public class DemoOrientBioAssembly {

    public static void main(String[] args){


        //String[] pdbIDs = new String[]{"4HHB","4AQ5","1LTI","1STP","4F88","2W6E","2LXC","3OE7","4INU","4D8s","4EAR","4IYQ","3ZKR"};

        String[] pdbIDs = new String[]{"4x2s"};

        int bioAssemblyNr = 1;

		/*
		    Local symmetry

			2WPD has 2 local symmetries.

			Other examples with a single local symmetry are:
			4F88 – local C8
			1LTI – local C5
			2W6E – local C3
			2LXC – local C2
			3OE7 – local C3

			Local Pseudosymmetry, structure only

			3ZDY
			3ZDX

			Helical

			1B47

		 */

        for ( String pdbID : pdbIDs)
        {
            try {

                runPDB(pdbID,bioAssemblyNr);

            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

    }

    public static void runPDB(String pdbID, int bioAssemblyNr) throws IOException, StructureException {



        pdbID = pdbID.toLowerCase();


        //Structure s = StructureIO.getBiologicalAssembly(pdbID, bioAssemblyNr);
        Structure s = readStructure(pdbID, bioAssemblyNr);

        QuatSymmetryParameters parameters = new QuatSymmetryParameters();
        parameters.setOnTheFly(true);
        parameters.setVerbose(true);


        CalcBioAssemblySymmetry calc = new CalcBioAssemblySymmetry(s, parameters);

        QuatSymmetryDetector detector = calc.orient();


        List<QuatSymmetryResults> globalResults = detector.getGlobalSymmetry();

        System.out.println("# of global results: " + globalResults.size());

        List<List<QuatSymmetryResults>> localResults = detector.getLocalSymmetries();



        showResults(s, pdbID + "[" + bioAssemblyNr + "] Global", globalResults);


        for (int counter = 0;counter < localResults.size() ; counter++){
            List<QuatSymmetryResults> localResultsL = localResults.get(counter);

            showResults(s,pdbID + "[" + bioAssemblyNr + "] Local #" + (counter+1) , localResultsL);
        }

        //determine default view:
        boolean defaultFound = false;

        for ( QuatSymmetryResults r : globalResults) {
            System.out.println(r.getSymmetry());
            //			if (! r.getRotationGroup().getPointGroup().equals("C1")) {
            //				defaultResult = r;
            //				defaultFound = true;
            //			} else if ( r.getSubunits().isPseudoSymmetric()) {
            //				defaultResult = r;
            //				defaultFound  = true;
            //			}
            if (r.isPreferredResult()) {

                defaultFound = true;
                System.out.println("!!!");
            }

        }
        if ( ! defaultFound) {
            for (List<QuatSymmetryResults> localResultSet : localResults) {
                for ( QuatSymmetryResults r : localResultSet) {
                    System.out.println(r.getSymmetry());
                    if ( r.isPreferredResult()) {

                    }
                }
            }
        }


    }

    private static void showResults(Structure s, String title,
                                    List<QuatSymmetryResults> results) {


        int count = 0 ;
        for (QuatSymmetryResults result: results) {

            String longTitle = title + " count:"+ count + " [" + result.getSubunits().getStoichiometry() +"]";

            String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; " +
                    "backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; " +
                    "color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;" ;
            count++;

            if ( result.getSubunits().isPseudoSymmetric()) {
                longTitle += " pseudosymmetric!";
            } else {
                System.out.println(" not pseudosymmetric!");

            }

            AxisAligner aligner = AxisAligner.getInstance(result);

            // use factory method to get point group specific instance of script generator
            JmolSymmetryScriptGenerator scriptGenerator = JmolSymmetryScriptGeneratorPointGroup.getInstance(aligner, "g");

            script += scriptGenerator.getOrientationWithZoom(0);
            script += scriptGenerator.drawPolyhedron();
            script += scriptGenerator.drawAxes();
            script += scriptGenerator.colorBySymmetry();


            longTitle += " M:" + result.getMethod();

            longTitle += String.format(" SEQ: %.2f - %.2f", result.getSubunits().getMinSequenceIdentity() ,result.getSubunits().getMaxSequenceIdentity());





            script += "draw axes* on; draw poly* on;";


            // show in Jmol...

            StructureAlignmentJmol jmol = new StructureAlignmentJmol();
            jmol.setStructure(s);

            jmol.setTitle(longTitle);
            jmol.evalString(script);
        }


    }



    private static Structure  readStructure(String pdbId, int bioAssemblyId) {
        // initialize the PDB_DIR env variable
        AtomCache cache = new AtomCache();
        cache.setUseMmCif(true);

        FileParsingParameters params = new FileParsingParameters();
        params.setStoreEmptySeqRes(true);
        params.setParseCAOnly(true);
        params.setLoadChemCompInfo(true);
        params.setAtomCaThreshold(Integer.MAX_VALUE);
        cache.setFileParsingParams(params);

        Structure structure = null;
        try {
            structure = StructureIO.getBiologicalAssembly(pdbId, bioAssemblyId);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (StructureException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        return structure;
    }
}
