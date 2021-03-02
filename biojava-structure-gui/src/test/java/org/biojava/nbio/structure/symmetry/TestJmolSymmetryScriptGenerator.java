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
package org.biojava.nbio.structure.symmetry;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.axis.RotationAxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.RotationGroup;
import org.biojava.nbio.structure.symmetry.core.Stoichiometry;
import org.biojava.nbio.structure.symmetry.core.SymmetryPerceptionMethod;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorDn;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author Spencer Bliven
 */
public class TestJmolSymmetryScriptGenerator {
    @Before
    public void setUp() {
    }

    @Test
    public void testPolygon() throws IOException, StructureException {
        Structure struc = StructureIO.getStructure("4hhb");
        QuatSymmetryParameters sp = new QuatSymmetryParameters();
		SubunitClustererParameters cp = new SubunitClustererParameters();

        QuatSymmetryResults results = QuatSymmetryDetector.calcGlobalSymmetry(struc, sp, cp);
        RotationAxisAligner axis = new RotationAxisAligner(results);
        JmolSymmetryScriptGeneratorDn gen = new JmolSymmetryScriptGeneratorDn(axis, "D3");

        String poly = gen.drawPolyhedron();
        String expected = "draw polyhedronD30 line{30.02,-39.95,0.59}{29.24,-0.53,40.00}{30.02,38.89,0.59}{30.80,-0.53,-38.82}{30.02,-39.95,0.59}{-30.00,-39.95,-0.60}{-30.79,-0.53,38.81}{-30.00,38.89,-0.60}{-29.22,-0.53,-40.01}{-30.00,-39.95,-0.60}width 0.45 color [x42ffd9] off;draw polyhedronD31 line{29.24,-0.53,40.00}{-30.79,-0.53,38.81}width 0.45 color [x42ffd9] off;draw polyhedronD32 line{30.02,38.89,0.59}{-30.00,38.89,-0.60}width 0.45 color [x42ffd9] off;draw polyhedronD33 line{30.80,-0.53,-38.82}{-29.22,-0.53,-40.01}width 0.45 color [x42ffd9] off;";
        assertEquals(expected, poly);
    }
}