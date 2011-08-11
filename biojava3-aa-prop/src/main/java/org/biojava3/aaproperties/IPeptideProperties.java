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
 * Created on 2011.05.09 by kohchuanhock
 *
 */
package org.biojava3.aaproperties;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.biojava3.aaproperties.xml.AminoAcidCompositionTable;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

/**
 *
 *
 *
 *
 *
 *
How To add a new executable

1) Add executable to the binaries folder. If it has source code and can be 
recompiled for different platforms include it under binaries/src 
Edit binaries/src setexecutableflag.sh and compilebin.sh scripts accordingly.

Added AAProperties.jar into jabaws/binaries/windows.
Need not edit *.sh since its a jar file.

2) Make sure all the dependencies of the software being installed are satisfied. 
If there are other binaries they should be included as well. Keep the dependant 
binaries in subfolder for the main executable. Update compile and setexecflag 
scripts accordingly.

Every dependencies are packed into the jar file hence nothing needs to be done here.

3) Make sure executable 
   - Does not have any hard links to its dependencies, e.g. is able to run from 
   any installation folder and does not contain any hard coded paths. 
   (TODO examples...)
   
Actually, there exists a hard link but maybe it is packaged into the jar file as well. That is the default XML Location. 
But I suppose it would be fine.

4) Describe executable in conf/Exectuable.properties. The lowercase name of the 
wrapper should be included in the name of the property for example Clustal 
properties all include clustal as a part of the name e.g. local.clustalw.bin
The same property for Mafft will be called local.mafft.bin. 

Added AAProperties description to conf/Executable.properties
Change AAProperties.jar to aaproperties.jar

GOOD TILL HERE

5) Add <ExecutableName>Limit.xml, <ExecutableName>Parameters.xml and 
<ExecutableName>Presets.xml. All are optional (should be at least). If the 
executable does not support parameters you do not have to refer to the 
XXXParameter.xml file into the Executable.properties file. The same is true for 
Presets and Limits. 


Added AAPropertiesLimits.xml and AAPropertiesParameters.xml
Ignore AAPropertiesPresets.xml for now

Question) Think that my AAPropertiesParameters.xml is not properly defined
i) My command prompt uses ' ' to separate instead of '='. Would that work?
ii) Some of my options do not have values. Would that be fine with the xml file?
iii) Any way to specify some options to be mandatory or optional?

6) Create a Java wrapper class for your executable. Create it within runner 
source directory. Examples of other wrappers can be found in compbio.runner.msa 
or compbio.runner.disorder packages. Wrapper should extend SkeletalExecutable<T> 
implements PipedExecutable<T> if you need to pass the input or collect the 
results from the standard in/out. Please see Mafft code as example. Wrapper 
should expend SkeletalExecutable<T> if input/output can be set as a parameter 
for an executable. Please see ClustalW code as example.

Created AAProperties.java under compbio.runner.sequence and extends SkeletalExecutable<AAProperties>

7) Create a testcase suit for your wrapper and run the test cases. 

Question) Stuck at Step 7.
i) What are the test cases are we supposed to run?
ii) Tried copying the JronnTester.testRunLocally but does not seems to work. Believe my step 6, Java wrapper class might not be working as well.

Initial attempt to provide AAProperties as a service in JABAWS

8) Create parser for the output files of your executable. Suggested location 
compbio.data.sequence.SequenceUtil  

9) Test the parser

Stop at step 9 then use the 

10) Decide which web services interface your executable is going to match. 
    For example if the executable output can be represented as 
    SequenceAnnotation then SequenceAnnotation interface might be appropriate. 
    For multiple sequence alignment an Msa interface should be used. 

11) If you find a web interface that matches your returning data type, then 
implement a web service which confirms to it within a webservices source folder 

12) Register web service in WEB-INF/ web.xml and sun-jaxws.xml

13) Add generated wsdl to wsbuild.xml ant script to generate the stubs

14) Run build-server task in wsbuild file. Watch for errors. If the task fails 
that means that JAXB cannot serialize some of the data structures. Add 
appropriate annotations to your data types.
Also check that 
  - you do not have interfaces to serialize. JAXB cannot serialize them.
  - you have a default no args constructor (can be private if you do not need it)
  - JAXB cannot serialize a Map, use custom data structure instead!
  - Enum cannot be serialized as its abstract class (do not confuse with enum 
  which is fine)
  - Fields serialization leave a little more space for manoeuvre, so use it. If 
  you do then you can accept and return interfaces, e.g. List, Map; abstract 
  classes etc, from your methods. 
  
If you have the data on the server side, but nothing is coming through to the 
client, this is a JAXB serialization problem. They tend to be very silent and 
thus hard to debug. Check your data structure can be serialized! 

13) Modify the client to work with your new web service. Update Services 
enumeration to include new service and ensure that all the methods of this 
enumeration take into account the new service. Update the client help text 
(client_help.txt) and insert it into the Constraints class.   

14) Test the web service with the client 

15) Test on the cluster...


 * 
 * TODO
 * 
 * Go into the websites to look at man_jaba_internals.html - Taught me how to link svn and abit about jaba internals.
 * ignore dependancies problem since i am using java
 * Jronn and AACon as examples
 * naming conventions are important
 * disorder.jronn
 * ignore clusters
 * goto SequenceUtil - Create readAAProperties
 * Be aware of the boundaries
 * 
 * runner.conservation.AACon.java
 * Download the WAR file and install the TomCat.
 * 
 * 
 * 
 * Look at Andreas email
 * http://www.biojava.org/download/maven/
 * http://emmy.rcsb.org:8080/cruisecontrol/
 * http://maven.apache.org/shared/maven-archiver/examples/classpath.html#Make
 * Commit only when you ensure it runs locally
 * Maven Compile around 20mins to 1hour after commit
 * 
 * distributionManagement
 * build
 * 
 * DONE
 * Change the return type of parseAAProp to ScoreManager
 * Adjust the configuration of pom.xml to generate a jar file with org.biojava3.aaproperties.CommandPrompt as the main class. However, need to rename it to AAProperties.jar
 * 
 * Question
 * Where to upload the jar file for the command prompt
 * 
 * An interface to generate some basic physico-chemical properties of protein sequences.<br/>
 * The following properties could be generated:
 * <p/>
 * Molecular weight<br/>
 * Absorbance<br/>
 * Extinction coefficient<br/>
 * Instability index<br/>
 * Apliphatic index<br/>
 * Average hydropathy value<br/>
 * Isoelectric point<br/>
 * Net charge at pH 7<br/>
 * Composition of specified amino acid<br/>
 * Composition of the 20 standard amino acid<br/>
 * @author kohchuanhock
 * @version 2011.05.09
 * @see PeptideProperties
 */
public interface IPeptideProperties{
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on <a href="http://web.expasy.org/findmod/findmod_masses.html">here</a>.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @see ProteinSequence
	 */
	public double getMolecularWeight(ProteinSequence sequence);
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the input files. These input files must be XML using the defined schema.
	 * Note that it assumes that ElementMass.xml file can be found in default location.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * 		xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 		xml file that details the composition of amino acids
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws JAXBException
	 * 		thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 		thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public double getMolecularWeight(ProteinSequence sequence, File aminoAcidCompositionFile) throws JAXBException, FileNotFoundException;
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the input files. These input files must be XML using the defined schema. 
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param elementMassFile
	 * 		xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 		xml file that details the composition of amino acids
	 * @return the total molecular weight of sequence + weight of water molecule
	 * @throws JAXBException
	 * 		thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 		thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public double getMolecularWeight(ProteinSequence sequence, File elementMassFile, File aminoAcidCompositionFile) 
		throws JAXBException, FileNotFoundException;
	
	/**
	 * Returns the molecular weight of sequence. The sequence argument must be a protein sequence consisting of only non-ambiguous characters.
	 * This method will sum the molecular weight of each amino acid in the
	 * sequence. Molecular weights are based on the AminoAcidCompositionTable. 
	 * Those input files must be XML using the defined schema.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCompositionTable
	 * 		a amino acid composition table obtained by calling IPeptideProperties.obtainAminoAcidCompositionTable
	 * @return the total molecular weight of sequence + weight of water molecule
	 */
	public double getMolecularWeightBasedOnXML(ProteinSequence sequence, AminoAcidCompositionTable aminoAcidCompositionTable);
	
	/**
	 * This method would initialize amino acid composition table based on the input xml files and stores the table for usage in future calls to 
	 * IPeptideProperties.getMolecularWeightBasedOnXML(ProteinSequence, AminoAcidCompositionTable).
	 * Note that ElementMass.xml is assumed to be able to be seen in default location.
	 * 
	 * @param aminoAcidCompositionFile
	 * 		xml file that details the composition of amino acids
	 * @return the initialized amino acid composition table
	 * @throws JAXBException
	 * 		thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 		thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File aminoAcidCompositionFile) 
		throws JAXBException, FileNotFoundException;
	
	/**
	 * This method would initialize amino acid composition table based on the input xml files and stores the table for usage in future calls to 
	 * IPeptideProperties.getMolecularWeightBasedOnXML(ProteinSequence, AminoAcidCompositionTable).
	 * 
	 * @param elementMassFile
	 * 		xml file that details the mass of each elements and isotopes
	 * @param aminoAcidCompositionFile
	 * 		xml file that details the composition of amino acids
	 * @return the initialized amino acid composition table
	 * @throws JAXBException
	 * 		thrown if unable to properly parse either elementMassFile or aminoAcidCompositionFile
	 * @throws FileNotFoundException
	 * 		thrown if either elementMassFile or aminoAcidCompositionFile are not found
	 */
	public AminoAcidCompositionTable obtainAminoAcidCompositionTable(File elementMassFile, File aminoAcidCompositionFile) 
		throws JAXBException, FileNotFoundException;

	/**
	 * Returns the extinction coefficient of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The extinction coefficient indicates how much light a protein absorbs at
	 * a certain wavelength. It is useful to have an estimation of this
	 * coefficient for following a protein which a spectrophotometer when
	 * purifying it. The computation of extinction coefficient follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 * 		true if Cys are assumed to be reduced and false if Cys are
	 *		assumed to form cystines
	 * @return the extinction coefficient of sequence
	 * @see ProteinSequence
	 */
	public double getExtinctionCoefficient(ProteinSequence sequence, boolean assumeCysReduced);

	/**
	 * Returns the absorbance (optical density) of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The computation of absorbance (optical density) follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param assumeCysReduced
	 * 		true if Cys are assumed to be reduced and false if Cys are
	 * 		assumed to form cystines
	 * @return the absorbance (optical density) of sequence
	 * @see ProteinSequence
	 */
	public double getAbsorbance(ProteinSequence sequence, boolean assumeCysReduced);
	
	/**
	 * Returns the instability index of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The instability index provides an estimate of the stability of your
	 * protein in a test tube. The computation of instability index follows the
	 * documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the instability index of sequence
	 * @see ProteinSequence
	 */
	public double getInstabilityIndex(ProteinSequence sequence);

	/**
	 * Returns the apliphatic index of sequence. The sequence argument must be a
	 * protein sequence consisting of only non-ambiguous characters.
	 * The aliphatic index of a protein is defined as the relative volume
	 * occupied by aliphatic side chains (alanine, valine, isoleucine, and
	 * leucine). It may be regarded as a positive factor for the increase of
	 * thermostability of globular proteins. The computation of aliphatic index
	 * follows the documentation in <a href="http://web.expasy.org/protparam/protparam-doc.html">here</a>.
	 * A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the aliphatic index of sequence
	 * @see ProteinSequence
	 */
	public double getApliphaticIndex(ProteinSequence sequence);

	/**
	 * Returns the average hydropathy value of sequence. The sequence argument
	 * must be a protein sequence consisting of only non-ambiguous characters.
	 * The average value for a sequence is calculated as the sum of hydropathy
	 * values of all the amino acids, divided by the number of residues in the
	 * sequence. Hydropathy values are based on (Kyte, J. and Doolittle, R.F.
	 * (1982) A simple method for displaying the hydropathic character of a
	 * protein. J. Mol. Biol. 157, 105-132).
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the average hydropathy value of sequence
	 * @see ProteinSequence
	 */
	public double getAvgHydropathy(ProteinSequence sequence);

	/**
	 * Returns the isoelectric point of sequence. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The isoelectric point is the pH at which the protein carries no net
	 * electrical charge. The isoelectric point will be computed based on
	 * approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#PI">here</a>
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the isoelectric point of sequence
	 * @see ProteinSequence
	 */
	public double getIsoelectricPoint(ProteinSequence sequence);

	/**
	 * Returns the net charge of sequence at pH 7. The sequence argument must be
	 * a protein sequence consisting of only non-ambiguous characters.
	 * The net charge will be computed using the approach stated in 
	 * <a href="http://www.innovagen.se/custom-peptide-synthesis/peptide-property-calculator/peptide-property-calculator-notes.asp#NetCharge>here</a>
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the net charge of sequence at pH 7
	 * @see ProteinSequence
	 */
	public double getNetCharge(ProteinSequence sequence);

	/**
	 * Returns the composition of specified amino acid in the sequence. The
	 * sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters. The aminoAcidCode must be a non-ambiguous
	 * character.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @param aminoAcidCode
	 * 		the code of the amino acid to compute
	 * @return the composition of specified amino acid in the sequence
	 * @see ProteinSequence 
	 * @see AminoAcidCompound
	 */
	public double getEnrichment(ProteinSequence sequence, AminoAcidCompound aminoAcidCode);

	/**
	 * Returns the composition of the 20 standard amino acid in the sequence.
	 * The sequence argument must be a protein sequence consisting of only
	 * non-ambiguous characters.
	 * The composition of an amino acid is the total number of its occurrence,
	 * divided by the total length of the sequence.
	 * 
	 * @param sequence
	 * 		a protein sequence consisting of non-ambiguous characters only
	 * @return the composition of the 20 standard amino acid in the sequence
	 * @see ProteinSequence 
	 * @see AminoAcidCompound
	 */
	public Map<AminoAcidCompound, Double> getAAComposition(ProteinSequence sequence);
}
