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
/**
 * <BODY>
 * <P>
 * Interfaces and classes for protein structure (PDB).
 * </P>
 * <p>
 * See also the <a href="https://github.com/biojava/biojava3-tutorial/blob/master/structure/README.md">BioJava 3 tutorial</a> for more information on the protein structure modules.
 * </p>
 * <h2>Parse PDB files</h2>
 * To load a PDB file see the <a href="io/PDBFileReader.html">PDBFileReader</a> class in the
 * <a href="io/package-summary.html">IO subpackage</a>.
 * <h2>Parse mmCif files</h2>
 * To laod a mmCif file see the <a href="io/MMCIFFileReader.html">MMCIFFileReader</a> class.
 * <h2>The Structure object</h2>
 * The <a href="Structure.html">Structure</a> object allows to access the PDB header information as well
 * as to the data from the ATOM records. The header information is currently available through the following objects:
 * <ul>
 * <li><a href="PDBHeader.html">PDBHeader</a></li>
 * <li><a href="DBRef.html">DBRef</a></li>
 * <li><a href="Compound.html">Compound</a></li>
 * </ul>
 * The structure object provides access to the data from the ATOM records through
 * a hierarchy of sub-object:
 * <pre>
 * <a href="Structure.html">Structure</a>
 * |
 * <a href="Chain.html">Chain</a>
 * |
 * <a href="Group.html">Group</a>
 * |
 * <a href="Atom.html">Atom</a>
 * </pre>
 * Learn more <a href="http://biojava.org/wiki/BioJava:CookBook:PDB:groups">how to work with groups</a>.
 * <h2>Other Features</h2>
 * <ul>
 * <li>Calculate <a href="align/ce/CeMain.html">protein structure alignments with CE and FATCAT.</a></li>
 * <li><a href="http://www.spice-3d.org/hibernatePDB/" target="_top">Serialize PDB files to databases using Hibernate</a></li>
 * <li><a href="Calc.html">Tools for performing calculations</a>
 * <li><a href="gui/BiojavaJmol.html">Display structures in Jmol</a></li>
 * </ul>
 * <p>
 * For more documentation on how to work with the Structure API please
 * see <a href="http://biojava.org/wiki/BioJava:CookBook#Protein_Structure" target="_top">
 * http://biojava.org/wiki/BioJava:CookBook#Protein_Structure</a>
 * </p>
 * </BODY>
 * @since 1.5
 */
package org.biojava.nbio.structure;